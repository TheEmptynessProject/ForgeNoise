const assert = require('assert');
const ForgeNoise = require('../forgeNoise.js');

describe('ForgeNoise', function() {
    const noise = new ForgeNoise(12345);

    it('should generate 2D noise within [-1, 1]', function() {
        for (let i = 0; i < 1000; i++) {
            const value = noise.generate2D(Math.random() * 100, Math.random() * 100);
            assert(value >= -1 && value <= 1, `Value ${value} is out of range`);
        }
    });

    it('should generate 3D noise within [-1, 1]', function() {
        for (let i = 0; i < 1000; i++) {
            const value = noise.generate3D(Math.random() * 100, Math.random() * 100, Math.random() * 100);
            assert(value >= -1 && value <= 1, `Value ${value} is out of range`);
        }
    });

    it('should generate consistent noise for the same input', function() {
        const value1 = noise.generate2D(1.5, 2.5);
        const value2 = noise.generate2D(1.5, 2.5);
        assert.strictEqual(value1, value2, 'Noise values should be consistent');
    });
});
