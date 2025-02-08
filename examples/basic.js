const ForgeNoise = require('../forgeNoise.js');

// Create a new noise generator with a seed
const noise = new ForgeNoise(12345);

// Generate some 2D noise
console.log('2D Noise:');
for (let y = 0; y < 5; y++) {
    for (let x = 0; x < 5; x++) {
        process.stdout.write(noise.generate2D(x * 0.1, y * 0.1).toFixed(2) + ' ');
    }
    console.log();
}

// Generate some 3D noise
console.log('\n3D Noise:');
for (let z = 0; z < 3; z++) {
    for (let y = 0; y < 3; y++) {
        for (let x = 0; x < 3; x++) {
            process.stdout.write(noise.generate3D(x * 0.1, y * 0.1, z * 0.1).toFixed(2) + ' ');
        }
        console.log();
    }
    console.log();
}
