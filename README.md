# ForgeNoise 🔥🗣

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

A high-performance JavaScript library for procedural noise generation, featuring multiple algorithms and fractal patterns. Perfect for games, visualizations, and creative coding projects.

## Examples

![Domain Warping](https://github.com/TheEmptynessProject/forgeNoise/blob/main/examples/DomainWarping.png)
![Ridged Multifractal](https://github.com/TheEmptynessProject/forgeNoise/blob/main/examples/RidgedMultifractal.png)
![Warped fBm with Turbulence](https://github.com/TheEmptynessProject/forgeNoise/blob/main/examples/WarpedfBmTurbulence.png)

## Features

- 🌪 Multiple noise algorithms: Perlin, Simplex, Worley, and Voronoi
- 🌀 Fractal noise generation (fBm, Ridged Multifractal)
- 🌐 Domain warping and turbulence effects
- 🧮 Seamless tiling patterns
- ⚡ Web-optimized performance
- 🌈 Customizable parameters for all noise types
- 🔢 Seedable randomness

## Installation

### CDN
```html
<script src="https://cdn.jsdelivr.net/gh/TheEmptynessProject/ForgeNoise@main/dist/forgeNoise.min.js"></script>
```

### Local
```html
<script src="path/to/forgeNoise.min.js"></script>
```

## Quick Start

```javascript
// Initialize generator
const noise = new ForgeNoise(seedNumber);

// Generate basic Perlin noise
const perlinValue = noise.generate2D(x, y);

// Create fractal noise
const fbmValue = noise.generateFractal2D(x, y, {
  octaves: 6,
  lacunarity: 2.0,
  persistence: 0.5
});

// Generate Worley cellular pattern
const worleyValue = noise.generateWorley2D01(x, y);

// Create domain-warped noise
const warped = noise.warp2D(x, y, {
  warpStrength: 2,
  warpScale: 1.5
});
```

## API Highlights

### Core Methods
- `new ForgeNoise([seed])` - Create new noise generator
- `generate2D(x, y)` - Classic Perlin noise (range: [-1, 1])
- `generateSimplex2D(x, y)` - Simplex noise implementation
- `generateWorley2D(x, y)` - Cellular/Worley noise

### Advanced Features
- `generateFractal2D()` - Fractional Brownian Motion
- `warp2D()` - Domain distortion effects
- `generateTiling2D()` - Seamless tiling patterns
- `generateVoronoi2D()` - Voronoi diagram generation

### Utility Methods
- `generate*01()` versions - Output mapped to [0, 1] range
- `setSeed(seed)` - Update generator seed
- `configure()` - Global noise parameters

[📚 Full API Documentation](https://yourapidocs.com)

## Examples

### Basic Perlin Noise
```javascript
const value = noise.generate2D01(x/20, y/20);
```

### Turbulent Terrain
```javascript
const height = noise.generateFractal2D(x/50, y/50, {
  octaves: 8,
  persistence: 0.65,
  turbulence: true
});
```

### Cellular Texturing
```javascript
const pattern = noise.generateWorley2D01(x/15, y/15) * 
               noise.generateSimplex2D01(x/30, y/30);
```

## Documentation

- [Interactive Examples](https://yourexamples.com)
- [API Reference](https://yourapidocs.com)

## Contributing

Contributions are welcome! Please read our 
[Contribution Guidelines](CONTRIBUTING.md) before submitting PRs.

## License

[MIT](LICENSE) © [TheEmptynessProject]
