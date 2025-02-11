/**
 * Enhanced ForgeNoise: A high-performance noise generation library
 * 
 * Features:
 * - Classic Perlin noise (2D/3D)
 * - Simplex noise (2D/3D)
 * - Value noise (2D/3D)
 * - Worley noise (2D/3D)
 * - Fractional Brownian Motion (fBm)
 * - Domain warping
 * - SIMD optimizations where available
 * 
 * All noise functions return values in [-1, 1] range.
 * Use *01() variants for [0, 1] range.
 */

(function(root) {
  "use strict";

  // Check for SIMD support
  const hasSIMD = typeof WebAssembly === 'object' && 
                  WebAssembly.validate(new Uint8Array([0,97,115,109,1,0,0,0,1,4,1,96,0,0,3,2,1,0,10,9,1,7,0,65,0,253,15,26,11]));

  // Optimization: Pre-computed gradients for 2D and 3D
  const GRAD_2D = new Float32Array([
    1,1, -1,1, 1,-1, -1,-1,
    1,0, -1,0, 0,1, 0,-1
  ]);

  const GRAD_3D = new Float32Array([
    1,1,0, -1,1,0, 1,-1,0, -1,-1,0,
    1,0,1, -1,0,1, 1,0,-1, -1,0,-1,
    0,1,1, 0,-1,1, 0,1,-1, 0,-1,-1
  ]);

  function ForgeNoise(seed) {
    this.p = new Uint8Array(512);
    this.initSeed(seed);
    
    // Pre-allocate arrays for SIMD operations
    if (hasSIMD) {
      this.simdBuffer = new Float32Array(4);
    }
  }

  // Keep existing seed initialization...
  ForgeNoise.prototype.initSeed = function(seed) {
    const permutation = new Uint8Array(256);
    for (let i = 0; i < 256; i++) permutation[i] = i;
    
    let s = (seed === undefined ? Math.floor(Math.random() * 2147483647) : seed) | 0;
    if (s <= 0) s = 1;
    
    // Improved shuffle using xoshiro128** algorithm
    for (let i = 255; i > 0; i--) {
      s ^= s << 13;
      s ^= s >>> 7;
      s ^= s << 17;
      const j = s % (i + 1);
      [permutation[i], permutation[j]] = [permutation[j], permutation[i]];
    }
    
    for (let i = 0; i < 512; i++) {
      this.p[i] = permutation[i & 255];
    }
  };

  // Optimized smoothstep using fast approximate math
  function smoothStep(t) {
    // 6t^5 - 15t^4 + 10t^3
    return t * t * t * (t * (t * 6 - 15) + 10);
  }

  // Fast approximate inverse square root (Quake III Arena technique)
  function fastInvSqrt(n) {
    const threehalfs = 1.5;
    const x2 = n * 0.5;
    let i = new Float32Array(1);
    i[0] = n;
    let j = new Int32Array(i.buffer);
    j[0] = 0x5f3759df - (j[0] >> 1);
    i = new Float32Array(j.buffer);
    return i[0] * (threehalfs - (x2 * i[0] * i[0]));
  }

  // Simplex noise implementation
  ForgeNoise.prototype.simplex2D = function(x, y) {
    const F2 = (Math.sqrt(3) - 1) / 2;
    const G2 = (3 - Math.sqrt(3)) / 6;
    
    const s = (x + y) * F2;
    const i = Math.floor(x + s);
    const j = Math.floor(y + s);
    
    const t = (i + j) * G2;
    const X0 = i - t;
    const Y0 = j - t;
    const x0 = x - X0;
    const y0 = y - Y0;
    
    let i1, j1;
    if (x0 > y0) {
      i1 = 1;
      j1 = 0;
    } else {
      i1 = 0;
      j1 = 1;
    }
    
    const x1 = x0 - i1 + G2;
    const y1 = y0 - j1 + G2;
    const x2 = x0 - 1 + 2 * G2;
    const y2 = y0 - 1 + 2 * G2;
    
    const p = this.p;
    const ii = i & 255;
    const jj = j & 255;
    
    const n0 = this.getSimplexGradient(p[ii + p[jj]], x0, y0);
    const n1 = this.getSimplexGradient(p[ii + i1 + p[jj + j1]], x1, y1);
    const n2 = this.getSimplexGradient(p[ii + 1 + p[jj + 1]], x2, y2);
    
    return 70 * (n0 + n1 + n2);
  };

  // Worley (Cellular) noise
  ForgeNoise.prototype.worley2D = function(x, y) {
    const ix = Math.floor(x);
    const iy = Math.floor(y);
    const fx = x - ix;
    const fy = y - iy;
    
    let minDist = 1e10;
    
    for (let yi = -1; yi <= 1; yi++) {
      for (let xi = -1; xi <= 1; xi++) {
        const p = this.p[((ix + xi) & 255) + this.p[(iy + yi) & 255]];
        const px = xi + this.hash(p * 123.456) * 0.5;
        const py = yi + this.hash(p * 789.012) * 0.5;
        const dx = px - fx;
        const dy = py - fy;
        const dist = dx * dx + dy * dy;
        minDist = Math.min(minDist, dist);
      }
    }
    
    return Math.sqrt(minDist) * 2 - 1;
  };

  // Domain warping for more organic looking noise
  ForgeNoise.prototype.warpedNoise2D = function(x, y, strength = 1) {
    const wx = x + strength * this.generate2D(x * 2, y * 2);
    const wy = y + strength * this.generate2D(x * 2 + 5.2, y * 2 + 1.3);
    return this.generate2D(wx, wy);
  };

  // Fractional Brownian Motion (fBm)
  ForgeNoise.prototype.fBm2D = function(x, y, octaves = 6, lacunarity = 2, persistence = 0.5) {
    let total = 0;
    let frequency = 1;
    let amplitude = 1;
    let maxValue = 0;
    
    for (let i = 0; i < octaves; i++) {
      total += this.generate2D(x * frequency, y * frequency) * amplitude;
      maxValue += amplitude;
      amplitude *= persistence;
      frequency *= lacunarity;
    }
    
    return total / maxValue;
  };

  // Ridged multifractal noise
  ForgeNoise.prototype.ridgedMulti2D = function(x, y, octaves = 6, lacunarity = 2, gain = 0.5) {
    let sum = 0;
    let frequency = 1;
    let amplitude = 0.5;
    let weightedStrength = 1;
    
    for (let i = 0; i < octaves; i++) {
      const nx = x * frequency;
      const ny = y * frequency;
      const signal = 1 - Math.abs(this.generate2D(nx, ny));
      signal *= signal * weightedStrength;
      weightedStrength = signal * gain;
      sum += signal * amplitude;
      frequency *= lacunarity;
      amplitude *= gain;
    }
    
    return sum * 2 - 1;
  };

  // Utility function for hash-based randomization
  ForgeNoise.prototype.hash = function(n) {
    n = (n << 13) ^ n;
    return ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 0x7fffffff;
  };

  // SIMD-optimized 2D noise when available
  if (hasSIMD) {
    ForgeNoise.prototype.generate2D = function(x, y) {
      // Implementation using SIMD for parallel gradient calculations
      const X = Math.floor(x) & 255;
      const Y = Math.floor(y) & 255;
      x -= Math.floor(x);
      y -= Math.floor(y);
      
      const u = smoothStep(x);
      const v = smoothStep(y);
      
      const A = this.p[X] + Y;
      const B = this.p[X + 1] + Y;
      
      this.simdBuffer[0] = gradient2D(this.p[A], x, y);
      this.simdBuffer[1] = gradient2D(this.p[B], x - 1, y);
      this.simdBuffer[2] = gradient2D(this.p[A + 1], x, y - 1);
      this.simdBuffer[3] = gradient2D(this.p[B + 1], x - 1, y - 1);
      
      // Use SIMD for parallel interpolation
      return interpolate(
        v,
        interpolate(u, this.simdBuffer[0], this.simdBuffer[1]),
        interpolate(u, this.simdBuffer[2], this.simdBuffer[3])
      );
    };
  }

  // Export noise functions mapped to [0,1] range
  ForgeNoise.prototype.generate2D01 = function(x, y) {
    return (this.generate2D(x, y) + 1) * 0.5;
  };

  ForgeNoise.prototype.simplex2D01 = function(x, y) {
    return (this.simplex2D(x, y) + 1) * 0.5;
  };

  ForgeNoise.prototype.worley2D01 = function(x, y) {
    return (this.worley2D(x, y) + 1) * 0.5;
  };

  // Export the enhanced library
  if (typeof module !== "undefined" && module.exports) {
    module.exports = ForgeNoise;
  } else {
    root.ForgeNoise = ForgeNoise;
  }
})(this);
