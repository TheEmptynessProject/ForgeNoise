# ForgeNoise

ForgeNoise is a high-performance Perlin noise generation library for JavaScript. It provides improved 2D and 3D noise generation with optimized algorithms.

## Installation

```
npm install forge-noise
```

## Usage

```
const ForgeNoise = require('forge-noise');

// Create a new noise generator (optionally with a seed)
const noise = new ForgeNoise(12345);

// Generate 2D noise
let value2D = noise.generate2D(3.14, 2.71);

// Generate 3D noise
let value3D = noise.generate3D(1.0, 2.0, 3.0);

// Generate noise in [0, 1] range
let value01 = noise.generate2D01(3.14, 2.71);
```

## API

- `new ForgeNoise(seed?)`: Create a new noise generator
- `generate2D(x, y)`: Generate 2D noise in [-1, 1] range
- `generate3D(x, y, z)`: Generate 3D noise in [-1, 1] range
- `generate2D01(x, y)`: Generate 2D noise in [0, 1] range
- `generate3D01(x, y, z)`: Generate 3D noise in [0, 1] range

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
