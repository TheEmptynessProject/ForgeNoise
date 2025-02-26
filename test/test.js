const assert = require('assert');
const ForgeNoise = require('../dist/forgeNoise.min.js');

describe('ForgeNoise', function() {
    this.timeout(5000);
    const noise = new ForgeNoise(12345);

    describe('Basic Noise Generation', function() {
        it('should generate 2D noise within [-1, 1]', function() {
            for (let i = 0; i < 1000; i++) {
                const value = noise.generate2D(Math.random() * 100, Math.random() * 100);
                assert(value >= -1 && value <= 1);
            }
        });

        it('should generate 3D noise within [-1, 1]', function() {
            for (let i = 0; i < 1000; i++) {
                const value = noise.generate3D(Math.random() * 100, Math.random() * 100, Math.random() * 100);
                assert(value >= -1 && value <= 1);
            }
        });

        it('should be consistent for same input', function() {
            const value1 = noise.generate2D(1.5, 2.5);
            const value2 = noise.generate2D(1.5, 2.5);
            assert.strictEqual(value1, value2);
        });

        it('should differ with different seeds', function() {
            const noise = new ForgeNoise(12345);
            const noise2 = new ForgeNoise(54321);
            const value1 = noise.generate2D(0.5, 0.5);
            const value2 = noise2.generate2D(0.5, 0.5);
            assert.notStrictEqual(value1, value2);
        });
    });

    describe('Simplex Noise', function() {
        it('should generate within [-1, 1]', function() {
            const noise = new ForgeNoise(12345);
            for (let i = 0; i < 100; i++) {
                const x = Math.random() * 100;
                const y = Math.random() * 100;
                const value = noise.generateSimplex2D(x, y);
                assert(value >= -1 && value <= 1);
            }
        });

        it('should be consistent', function() {
            const noise = new ForgeNoise(12345);
            const x = 1.5, y = 2.5;
            const value1 = noise.generateSimplex2D(x, y);
            const value2 = noise.generateSimplex2D(x, y);
            assert.strictEqual(value1, value2);
        });
    });

    describe('Worley Noise', function() {
        it('should generate within [-1, 1]', function() {
            const noise = new ForgeNoise(12345);
            for (let i = 0; i < 100; i++) {
                const x = Math.random() * 100;
                const y = Math.random() * 100;
                const value = noise.generateWorley2D(x, y) * 2 - 1;
                assert(value >= -1 && value <= 1);
            }
        });

        it('should be consistent', function() {
            const noise = new ForgeNoise(12345);
            const x = 1.5, y = 2.5;
            const value1 = noise.generateWorley2D(x, y);
            const value2 = noise.generateWorley2D(x, y);
            assert.strictEqual(value1, value2);
        });
    });

    describe('Domain Warping', function() {
        it('should generate within [-1, 1]', function() {
            const noise = new ForgeNoise(12345);
            for (let i = 0; i < 100; i++) {
                const x = Math.random() * 100;
                const y = Math.random() * 100;
                const value = noise.warp2D(x, y, { warpStrength: 1.0 });
                assert(value >= -1 && value <= 1);
            }
        });

        it('should respect strength', function() { //Precision point errors
            const noise = new ForgeNoise(12345);
            const x = 1.5, y = 2.5;
            const baseValue = noise.warp2D(x, y, { warpStrength: 0, warpScale: 10 });
            const warpedValue = noise.warp2D(x, y, { warpStrength: 1, warpScale: 10 });
            assert.strictEqual(baseValue, noise.generate2D(x, y));
            assert.notStrictEqual(baseValue, warpedValue);
        });
    });

    describe('Fractional Brownian Motion', function() {
        it('should generate within [-1, 1]', function() {
            const noise = new ForgeNoise(12345);
            for (let i = 0; i < 100; i++) {
                const x = Math.random() * 100;
                const y = Math.random() * 100;
                const value = noise.generateFractal2D(x, y, { octaves: 4 });
                assert(value >= -1 && value <= 1);
            }
        });

        it('should respect octaves', function() {
            const noise = new ForgeNoise(12345);
            const x = 0.5, y = 0.5;
            const value1 = noise.generateFractal2D(x, y, { octaves: 1 });
            const value2 = noise.generateFractal2D(x, y, { octaves: 4 });
            assert.notStrictEqual(value1, value2);
        });
    });

    describe('Ridged Multifractal', function() {
        it('should generate within [-1, 1]', function() {
            const noise = new ForgeNoise(12345);
            for (let i = 0; i < 100; i++) {
                const x = Math.random() * 100;
                const y = Math.random() * 100;
                const value = noise.generateFractal2D(x, y, { octaves: 4, ridged: true });
                assert(value >= -1 && value <= 1);
            }
        });

        it('should respect gain', function() {
            const noise = new ForgeNoise(12345);
            const x = 0.5, y = 0.5;
            const value1 = noise.generateFractal2D(x, y, { octaves: 4, ridged: true, persistence: 0.5 });
            const value2 = noise.generateFractal2D(x, y, { octaves: 4, ridged: true, persistence: 0.8 });
            assert.notStrictEqual(value1, value2);
        });
    });

    describe('Range Mapping', function() {
        it('should map 2D to [0, 1]', function() {
            for (let i = 0; i < 1000; i++) {
                const value = noise.generate2D01(Math.random() * 100, Math.random() * 100);
                assert(value >= 0 && value <= 1);
            }
        });

        it('should map Simplex to [0, 1]', function() {
            const noise = new ForgeNoise(12345);
            for (let i = 0; i < 100; i++) {
                const x = Math.random() * 100;
                const y = Math.random() * 100;
                const value = noise.generateSimplex2D01(x, y);
                assert(value >= 0 && value <= 1);
            }
        });

        it('should map Worley to [0, 1]', function() {
            const noise = new ForgeNoise(12345);
            for (let i = 0; i < 100; i++) {
                const x = Math.random() * 100;
                const y = Math.random() * 100;
                const value = noise.generateWorley2D01(x, y);
                assert(value >= 0 && value <= 1);
            }
        });
    });

    describe('Performance', function() {
        it('should generate 2D noise efficiently', function() {
            const start = performance.now();
            for (let i = 0; i < 10000; i++) noise.generate2D(Math.random() * 100, Math.random() * 100);
            const time = (performance.now() - start) / 10000;
            console.log(`2D noise: ${time.toFixed(3)}ms per call`);
            assert(time < 1);
        });

        it('should generate Simplex noise efficiently', function() {
            const noise = new ForgeNoise(12345);
            const start = performance.now();
            for (let i = 0; i < 10000; i++) {
                noise.generateSimplex2D(i * 0.1, i * 0.2);
            }
            const time = performance.now() - start;
            assert(time < 50);
        });
    });

    describe('Edge Cases', function() {
        it('should handle integers', function() {
            const value = noise.generate2D(1, 1);
            assert(typeof value === 'number' && !isNaN(value));
        });

        it('should handle negatives', function() {
            const value = noise.generate2D(-1.5, -2.5);
            assert(value >= -1 && value <= 1);
        });

        it('should handle zeros', function() {
            const value = noise.generate2D(0, 0);
            assert(typeof value === 'number' && !isNaN(value));
        });

        it('should handle large coordinates', function() {
            const value = noise.generate2D(1e6, 1e6);
            assert(value >= -1 && value <= 1);
        });
    });
});
