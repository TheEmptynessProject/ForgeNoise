document.addEventListener('DOMContentLoaded', function() {

    window.ForgeNoise = ForgeNoise;

    window.noise = new ForgeNoise(12345); 

    const canvases = document.querySelectorAll('canvas');
    const preBlocks = document.querySelectorAll('pre[contenteditable="true"]');
    const seedBlock = document.getElementById('noise-seed');
    let zoomedCanvas = null;
    let zoomLevel = 1;
    let isDragging = false;
    let startX, startY, translateX = 0, translateY = 0;
    let dragStartX, dragStartY, movedDistance = 0;

    const overlay = document.createElement('div');
    overlay.className = 'zoom-overlay';
    overlay.style.position = 'fixed';
    overlay.style.top = '0';
    overlay.style.left = '0';
    overlay.style.width = '100%';
    overlay.style.height = '100%';
    overlay.style.background = 'rgba(0, 0, 0, 0.5)';
    overlay.style.zIndex = '999';
    overlay.style.display = 'none';
    document.body.appendChild(overlay);

    function renderCanvas(canvas, code) {
        const ctx = canvas.getContext('2d');
        const width = canvas.width;
        const height = canvas.height;
        const imageData = ctx.createImageData(width, height);

        try {
            if (!window.noise) throw new Error('Noise instance is not defined');
            const getValue = new Function('x', 'y', 'noise', `
                return (function() {
                    ${code}
                })();
            `);
            getValue(0, 0, window.noise);
            for (let x = 0; x < width; x++) {
                for (let y = 0; y < height; y++) {
                    const value = getValue(x, y, window.noise);
                    const index = (x + y * width) * 4;
                    const color = Math.floor(Math.min(Math.max(value, 0), 1) * 255);
                    imageData.data[index] = color;
                    imageData.data[index + 1] = color;
                    imageData.data[index + 2] = color;
                    imageData.data[index + 3] = 255;
                }
            }
            ctx.putImageData(imageData, 0, 0);
        } catch (e) {
            console.error(`Error rendering ${canvas.id}: ${e.message}`);
            ctx.fillStyle = '#ffcccc';
            ctx.fillRect(0, 0, width, height);
            ctx.fillStyle = '#ff0000';
            ctx.font = '14px Arial';
            ctx.fillText(`Error: ${e.message}`, 10, 20);
        }
    }

    function updateAllCanvases() {
        preBlocks.forEach(pre => {
            if (pre.id !== 'noise-seed') {
                const canvasId = pre.getAttribute('data-canvas');
                const canvas = document.getElementById(canvasId);
                const code = pre.textContent.trim();
                renderCanvas(canvas, code);
            }
        });
    }

    seedBlock.addEventListener('input', () => {
        const seedCode = seedBlock.textContent.trim();
        window.noise = undefined;
        try {
            if (seedCode) {
                eval(seedCode);
            }
            if (!window.noise) throw new Error('Noise instance is not defined');
            seedBlock.style.backgroundColor = '';
            updateAllCanvases();
        } catch (e) {
            console.error('Error updating noise seed:', e);
            seedBlock.style.backgroundColor = '#ffcccc';
            updateAllCanvases();
        }
    });

    eval(seedBlock.textContent.trim());
    updateAllCanvases();

    preBlocks.forEach(pre => {
        if (pre.id !== 'noise-seed') {
            pre.addEventListener('input', () => {
                const canvasId = pre.getAttribute('data-canvas');
                const canvas = document.getElementById(canvasId);
                const code = pre.textContent.trim();
                renderCanvas(canvas, code);
            });
        }
    });

    canvases.forEach(canvas => {

        canvas.style.cursor = 'zoom-in';
        canvas.addEventListener('mousemove', (e) => {
            canvas.style.cursor = e.ctrlKey ? 'zoom-out' : 'zoom-in';
        });
        canvas.addEventListener('mousedown', () => {
            canvas.style.cursor = window.event.ctrlKey ? 'zoom-out' : 'zoom-in';
        });
        canvas.addEventListener('mouseup', () => {
            canvas.style.cursor = 'zoom-in';
        });

        canvas.addEventListener('click', (e) => {
            e.stopPropagation();
            if (movedDistance > 5) return; 
            const ctrlPressed = e.ctrlKey;
            if (zoomedCanvas === canvas) {

                zoomLevel = ctrlPressed ? Math.max(1, zoomLevel - 0.5) : zoomLevel + 0.5;
                canvas.style.transform = `translate(${translateX}px, ${translateY}px) scale(${zoomLevel})`;
            } else {
                if (zoomedCanvas) {
                    revertCanvas();
                }

                canvas.style.position = 'fixed';
                canvas.style.width = '1024px';
                canvas.style.height = '1024px';
                canvas.style.top = '50%';
                canvas.style.left = '50%';
                canvas.style.transform = 'translate(-50%, -50%)';
                canvas.style.zIndex = '1000';
                zoomLevel = 1;
                translateX = 0;
                translateY = 0;
                zoomedCanvas = canvas;

                overlay.style.display = 'block';
                document.body.style.overflow = 'hidden';

                let closeBtn = canvas.nextElementSibling;
                if (!closeBtn || !closeBtn.classList.contains('close-btn')) {
                    closeBtn = document.createElement('button');
                    closeBtn.textContent = 'Close';
                    closeBtn.className = 'close-btn';
                    closeBtn.style.position = 'fixed'; 
                    closeBtn.style.top = '10px';
                    closeBtn.style.right = '10px';
                    closeBtn.style.background = '#e74c3c';
                    closeBtn.style.color = 'white';
                    closeBtn.style.border = 'none';
                    closeBtn.style.padding = '5px 10px';
                    closeBtn.style.borderRadius = '4px';
                    closeBtn.style.cursor = 'pointer';
                    closeBtn.style.zIndex = '1001'; 
                    canvas.parentNode.insertBefore(closeBtn, canvas.nextSibling);
                    closeBtn.addEventListener('click', (e) => {
                        e.stopPropagation();
                        revertCanvas();
                    });
                }
                closeBtn.style.display = 'block';
            }
        });

        canvas.addEventListener('dblclick', (e) => {
            e.stopPropagation();
            if (zoomedCanvas === canvas) {
                revertCanvas();
            }
        });

        canvas.addEventListener('mousedown', (e) => {
            if (zoomedCanvas === canvas) {
                isDragging = true;
                dragStartX = e.clientX;
                dragStartY = e.clientY;
                startX = e.clientX - translateX;
                startY = e.clientY - translateY;
                canvas.style.cursor = 'grab';
            }
        });

        canvas.addEventListener('mousemove', (e) => {
            if (isDragging && zoomedCanvas === canvas) {
                translateX = e.clientX - startX;
                translateY = e.clientY - startY;
                movedDistance = Math.sqrt(
                    Math.pow(e.clientX - dragStartX, 2) + Math.pow(e.clientY - dragStartY, 2)
                );
                canvas.style.transform = `translate(${translateX}px, ${translateY}px) scale(${zoomLevel})`;
            }
        });

        canvas.addEventListener('mouseup', () => {
            if (isDragging && zoomedCanvas === canvas) {
                isDragging = false;
                canvas.style.cursor = 'zoom-in';
                setTimeout(() => { movedDistance = 0; }, 100);
            }
        });

        canvas.addEventListener('mouseleave', () => {
            if (isDragging && zoomedCanvas === canvas) {
                isDragging = false;
                canvas.style.cursor = 'zoom-in';
                movedDistance = 0;
            }
        });
    });

    function revertCanvas() {
        if (zoomedCanvas) {
            zoomedCanvas.style.position = '';
            zoomedCanvas.style.width = '256px';
            zoomedCanvas.style.height = '256px';
            zoomedCanvas.style.top = '';
            zoomedCanvas.style.left = '';
            zoomedCanvas.style.transform = '';
            zoomedCanvas.style.zIndex = '';
            zoomedCanvas.style.cursor = 'zoom-in';
            zoomedCanvas.nextElementSibling.style.display = 'none';
            zoomedCanvas = null;
            zoomLevel = 1;
            translateX = 0;
            translateY = 0;
            overlay.style.display = 'none';
            document.body.style.overflow = '';
        }
    }

    document.addEventListener('click', (e) => {
        if (zoomedCanvas && !zoomedCanvas.contains(e.target) && e.target !== seedBlock && !e.target.classList.contains('close-btn')) {
            revertCanvas();
        }
    });

    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && zoomedCanvas) {
            revertCanvas();
        }
    });
});