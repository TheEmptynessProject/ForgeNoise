/**
 * ForgeNoise: A custom, high‑performance noise library.
 *
 * Provides improved Perlin noise in 2D and 3D.
 * Noise values are in the range [–1, +1]. (Use noise2D01()/noise3D01() to remap to [0, 1].)
 *
 * Usage:
 *   // With a seed (optional)
 *   const noise = new ForgeNoise(12345);
 *   let value2D = noise.generate2D(3.14, 2.71);
 *   let value01 = noise.generate2D01(3.14, 2.71);
 *
 *   // In a module system:
 *   // module.exports = ForgeNoise;
 */

(function(root) {
  "use strict";

  // The constructor accepts an optional numeric seed.
  function ForgeNoise(seed) {
    // p: permutation table (512 values) for fast indexing.
    this.p = new Uint8Array(512);
    this.initSeed(seed);
  }

  /**
   * Initialize the noise generator with a seed.
   * A simple linear congruential generator (LCG) shuffles the initial permutation.
   * If no seed is provided, Math.random() is used.
   */
  ForgeNoise.prototype.initSeed = function(seed) {
    // Create a 256-element permutation array with values 0..255.
    const permutation = new Uint8Array(256);
    for (let i = 0; i < 256; i++) {
      permutation[i] = i;
    }
    // Use the seed (or a random number) to shuffle.
    // Here we use the minimal standard LCG: s = (s * 16807) % 2147483647.
    let s = (seed === undefined ? Math.floor(Math.random() * 2147483647) : seed) | 0;
    if (s <= 0) s = 1; // avoid zero
    for (let i = 255; i > 0; i--) {
      s = (s * 16807) % 2147483647;
      // Get a pseudo‑random index j such that 0 ≤ j ≤ i.
      const j = s % (i + 1);
      // Swap permutation[i] and permutation[j].
      const tmp = permutation[i];
      permutation[i] = permutation[j];
      permutation[j] = tmp;
    }
    // Duplicate the permutation array into this.p (of length 512)
    // so that we avoid overflow when using (index & 255).
    for (let i = 0; i < 512; i++) {
      this.p[i] = permutation[i & 255];
    }
  };

  // Fade (smoothing) function: 6t^5 – 15t^4 + 10t^3.
  // This quintic polynomial eases coordinate values so that they ease towards integral values.
  function smoothStep(t) {
    return t * t * t * (t * (t * 6 - 15) + 10);
  }

  // Linear interpolation between a and b by fraction t.
  function interpolate(t, a, b) {
    return a + t * (b - a);
  }

  // 2D gradient function.
  // The lower 3 bits of the hash determine one of 8 possible directions.
  function gradient2D(hash, x, y) {
    const h = hash & 7; // 0..7
    // Choose x or y as the first component based on h.
    const u = h < 4 ? x : y;
    const v = h < 4 ? y : x;
    // Return the dot product with a sign change based on the lower bits.
    return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
  }

  // 3D gradient function.
  // Uses the lower 4 bits to choose from 16 possible directions.
  function gradient3D(hash, x, y, z) {
    const h = hash & 15; // 0..15
    const u = h < 8 ? x : y;
    const v = h < 4 ? y : (h === 12 || h === 14 ? x : z);
    return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
  }

  /**
   * generate2D(x, y): Returns a noise value for the point (x, y) in 2D.
   * The output is in the range [–1, +1].
   */
  ForgeNoise.prototype.generate2D = function(x, y) {
    // Find unit grid cell containing point.
    const X = Math.floor(x) & 255;
    const Y = Math.floor(y) & 255;
    // Find relative x, y in cell.
    x = x - Math.floor(x);
    y = y - Math.floor(y);
    // Compute fade curves for x and y.
    const u = smoothStep(x);
    const v = smoothStep(y);
    const p = this.p;
    const A = p[X] + Y;
    const B = p[X + 1] + Y;
    // Interpolate the results from the four corners of the cell.
    const result = interpolate(
      v,
      interpolate(u, gradient2D(p[A], x, y), gradient2D(p[B], x - 1, y)),
      interpolate(u, gradient2D(p[A + 1], x, y - 1), gradient2D(p[B + 1], x - 1, y - 1))
    );
    return result;
  };

  /**
   * generate3D(x, y, z): Returns a noise value for the point (x, y, z) in 3D.
   * The output is in the range [–1, +1].
   */
  ForgeNoise.prototype.generate3D = function(x, y, z) {
    const X = Math.floor(x) & 255;
    const Y = Math.floor(y) & 255;
    const Z = Math.floor(z) & 255;
    x = x - Math.floor(x);
    y = y - Math.floor(y);
    z = z - Math.floor(z);
    const u = smoothStep(x);
    const v = smoothStep(y);
    const w = smoothStep(z);
    const p = this.p;
    const A  = p[X] + Y,  AA = p[A] + Z, AB = p[A + 1] + Z;
    const B  = p[X + 1] + Y, BA = p[B] + Z, BB = p[B + 1] + Z;
    const result = interpolate(
      w,
      interpolate(
        v,
        interpolate(u, gradient3D(p[AA], x, y, z), gradient3D(p[BA], x - 1, y, z)),
        interpolate(u, gradient3D(p[AB], x, y - 1, z), gradient3D(p[BB], x - 1, y - 1, z))
      ),
      interpolate(
        v,
        interpolate(u, gradient3D(p[AA + 1], x, y, z - 1), gradient3D(p[BA + 1], x - 1, y, z - 1)),
        interpolate(u, gradient3D(p[AB + 1], x, y - 1, z - 1), gradient3D(p[BB + 1], x - 1, y - 1, z - 1))
      )
    );
    return result;
  };

  /**
   * Convenience helpers to map noise output from [–1, +1] to [0, 1].
   */
  ForgeNoise.prototype.generate2D01 = function(x, y) {
    return (this.generate2D(x, y) + 1) * 0.5;
  };

  ForgeNoise.prototype.generate3D01 = function(x, y, z) {
    return (this.generate3D(x, y, z) + 1) * 0.5;
  };

  // Export the library for CommonJS, AMD, or as a global.
  if (typeof module !== "undefined" && module.exports) {
    module.exports = ForgeNoise;
  } else {
    root.ForgeNoise = ForgeNoise;
  }
})(this);
