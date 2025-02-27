(function(root) {
	"use strict";

	const grad3 = [
		[1, 1, 0],
		[-1, 1, 0],
		[1, -1, 0],
		[-1, -1, 0],
		[1, 0, 1],
		[-1, 0, 1],
		[1, 0, -1],
		[-1, 0, -1],
		[0, 1, 1],
		[0, -1, 1],
		[0, 1, -1],
		[0, -1, -1]
	];
	const grad2 = grad3.map(v => [v[0], v[1]]);

	const grad4 = [
		[0, 1, 1, 1],
		[0, 1, 1, -1],
		[0, 1, -1, 1],
		[0, 1, -1, -1],
		[0, -1, 1, 1],
		[0, -1, 1, -1],
		[0, -1, -1, 1],
		[0, -1, -1, -1],
		[1, 0, 1, 1],
		[1, 0, 1, -1],
		[1, 0, -1, 1],
		[1, 0, -1, -1],
		[-1, 0, 1, 1],
		[-1, 0, 1, -1],
		[-1, 0, -1, 1],
		[-1, 0, -1, -1],
		[1, 1, 0, 1],
		[1, 1, 0, -1],
		[1, -1, 0, 1],
		[1, -1, 0, -1],
		[-1, 1, 0, 1],
		[-1, 1, 0, -1],
		[-1, -1, 0, 1],
		[-1, -1, 0, -1],
		[1, 1, 1, 0],
		[1, 1, -1, 0],
		[1, -1, 1, 0],
		[1, -1, -1, 0],
		[-1, 1, 1, 0],
		[-1, 1, -1, 0],
		[-1, -1, 1, 0],
		[-1, -1, -1, 0]
	];

	function dot2D(g, x, y) {
		return g[0] * x + g[1] * y;
	}

	function dot3D(g, x, y, z) {
		return g[0] * x + g[1] * y + g[2] * z;
	}

	function dot4D(g, x, y, z, w) {
		return g[0] * x + g[1] * y + g[2] * z + g[3] * w;
	}

	function smoothStep(t) {
		return t * t * t * (t * (t * 6 - 15) + 10);
	}

	function interpolate(t, a, b) {
		return a + t * (b - a);
	}

	function ForgeNoise(seed) {
		this.p = new Uint8Array(512);
		this.initSeed(seed);
	}

	ForgeNoise.prototype.initSeed = function(seed) {
		const permutation = new Uint8Array(256);
		for (let i = 0; i < 256; i++) {
			permutation[i] = i;
		}
		let s = (seed === undefined ? Math.floor(Math.random() * 2147483647) : seed) | 0;
		if (s <= 0) s = 1;
		for (let i = 255; i > 0; i--) {
			s = (s * 16807) % 2147483647;
			const j = s % (i + 1);
			[permutation[i], permutation[j]] = [permutation[j], permutation[i]];
		}
		for (let i = 0; i < 512; i++) {
			this.p[i] = permutation[i & 255];
		}
	};

	function gradient2D(hash, x, y) {
		const h = hash & 7;
		const u = h < 4 ? x : y;
		const v = h < 4 ? y : x;
		return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
	}

	function gradient3D(hash, x, y, z) {
		const h = hash & 15;
		const u = h < 8 ? x : y;
		const v = h < 4 ? y : (h === 12 || h === 14 ? x : z);
		return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
	}

	ForgeNoise.prototype.generate2D = function(x, y) {
		const X = Math.floor(x) & 255;
		const Y = Math.floor(y) & 255;
		x -= Math.floor(x);
		y -= Math.floor(y);
		const u = smoothStep(x);
		const v = smoothStep(y);
		const p = this.p;
		const A = p[X] + Y;
		const B = p[X + 1] + Y;
		return interpolate(
			v,
			interpolate(u, gradient2D(p[A], x, y), gradient2D(p[B], x - 1, y)),
			interpolate(u, gradient2D(p[A + 1], x, y - 1), gradient2D(p[B + 1], x - 1, y - 1))
		);
	};

	ForgeNoise.prototype.generate3D = function(x, y, z) {
		const X = Math.floor(x) & 255;
		const Y = Math.floor(y) & 255;
		const Z = Math.floor(z) & 255;
		x -= Math.floor(x);
		y -= Math.floor(y);
		z -= Math.floor(z);
		const u = smoothStep(x);
		const v = smoothStep(y);
		const w = smoothStep(z);
		const p = this.p;
		const A = p[X] + Y,
			AA = p[A] + Z,
			AB = p[A + 1] + Z;
		const B = p[X + 1] + Y,
			BA = p[B] + Z,
			BB = p[B + 1] + Z;
		return interpolate(
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
	};

	ForgeNoise.prototype.generateSimplex2D = function(x, y) {
		const F2 = 0.5 * (Math.sqrt(3) - 1);
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
		const ii = i & 255;
		const jj = j & 255;
		const p = this.p;
		const gi0 = p[ii + p[jj]] % 12;
		const gi1 = p[ii + i1 + p[jj + j1]] % 12;
		const gi2 = p[ii + 1 + p[jj + 1]] % 12;
		let n0 = 0,
			n1 = 0,
			n2 = 0;
		let t0 = 0.5 - x0 * x0 - y0 * y0;
		if (t0 >= 0) {
			t0 *= t0;
			n0 = t0 * t0 * dot2D(grad2[gi0], x0, y0);
		}
		let t1 = 0.5 - x1 * x1 - y1 * y1;
		if (t1 >= 0) {
			t1 *= t1;
			n1 = t1 * t1 * dot2D(grad2[gi1], x1, y1);
		}
		let t2 = 0.5 - x2 * x2 - y2 * y2;
		if (t2 >= 0) {
			t2 *= t2;
			n2 = t2 * t2 * dot2D(grad2[gi2], x2, y2);
		}
		return 70 * (n0 + n1 + n2);
	};

	ForgeNoise.prototype.generateSimplex3D = function(x, y, z) {
		const F3 = 1 / 3;
		const G3 = 1 / 6;
		const s = (x + y + z) * F3;
		const i = Math.floor(x + s);
		const j = Math.floor(y + s);
		const k = Math.floor(z + s);
		const t = (i + j + k) * G3;
		const X0 = i - t;
		const Y0 = j - t;
		const Z0 = k - t;
		const x0 = x - X0;
		const y0 = y - Y0;
		const z0 = z - Z0;
		let i1, j1, k1, i2, j2, k2;
		if (x0 >= y0) {
			if (y0 >= z0) {
				i1 = 1;
				j1 = 0;
				k1 = 0;
				i2 = 1;
				j2 = 1;
				k2 = 0;
			} else if (x0 >= z0) {
				i1 = 1;
				j1 = 0;
				k1 = 0;
				i2 = 1;
				j2 = 0;
				k2 = 1;
			} else {
				i1 = 0;
				j1 = 0;
				k1 = 1;
				i2 = 1;
				j2 = 0;
				k2 = 1;
			}
		} else {
			if (y0 < z0) {
				i1 = 0;
				j1 = 0;
				k1 = 1;
				i2 = 0;
				j2 = 1;
				k2 = 1;
			} else if (x0 < z0) {
				i1 = 0;
				j1 = 1;
				k1 = 0;
				i2 = 0;
				j2 = 1;
				k2 = 1;
			} else {
				i1 = 0;
				j1 = 1;
				k1 = 0;
				i2 = 1;
				j2 = 1;
				k2 = 0;
			}
		}
		const x1 = x0 - i1 + G3;
		const y1 = y0 - j1 + G3;
		const z1 = z0 - k1 + G3;
		const x2 = x0 - i2 + 2 * G3;
		const y2 = y0 - j2 + 2 * G3;
		const z2 = z0 - k2 + 2 * G3;
		const x3 = x0 - 1 + 3 * G3;
		const y3 = y0 - 1 + 3 * G3;
		const z3 = z0 - 1 + 3 * G3;
		const ii = i & 255;
		const jj = j & 255;
		const kk = k & 255;
		const p = this.p;
		const gi0 = p[ii + p[jj + p[kk]]] % 12;
		const gi1 = p[ii + i1 + p[jj + j1 + p[kk + k1]]] % 12;
		const gi2 = p[ii + i2 + p[jj + j2 + p[kk + k2]]] % 12;
		const gi3 = p[ii + 1 + p[jj + 1 + p[kk + 1]]] % 12;
		let n0 = 0,
			n1 = 0,
			n2 = 0,
			n3 = 0;
		let t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
		if (t0 >= 0) {
			t0 *= t0;
			n0 = t0 * t0 * dot3D(grad3[gi0], x0, y0, z0);
		}
		let t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
		if (t1 >= 0) {
			t1 *= t1;
			n1 = t1 * t1 * dot3D(grad3[gi1], x1, y1, z1);
		}
		let t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
		if (t2 >= 0) {
			t2 *= t2;
			n2 = t2 * t2 * dot3D(grad3[gi2], x2, y2, z2);
		}
		let t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
		if (t3 >= 0) {
			t3 *= t3;
			n3 = t3 * t3 * dot3D(grad3[gi3], x3, y3, z3);
		}
		return 32 * (n0 + n1 + n2 + n3);
	};

	ForgeNoise.prototype.generateValue2D = function(x, y) {
		const i0 = Math.floor(x);
		const i1 = i0 + 1;
		const j0 = Math.floor(y);
		const j1 = j0 + 1;
		const fx = x - i0;
		const fy = y - j0;
		const u = smoothStep(fx);
		const v = smoothStep(fy);
		const p = this.p;
		const a00 = p[(i0 + p[j0 & 255]) & 255] / 255;
		const a01 = p[(i0 + p[j1 & 255]) & 255] / 255;
		const a10 = p[(i1 + p[j0 & 255]) & 255] / 255;
		const a11 = p[(i1 + p[j1 & 255]) & 255] / 255;
		const mixX0 = interpolate(u, a00, a10);
		const mixX1 = interpolate(u, a01, a11);
		return interpolate(v, mixX0, mixX1);
	};

	ForgeNoise.prototype.generateValue3D = function(x, y, z) {
		const i0 = Math.floor(x);
		const i1 = i0 + 1;
		const j0 = Math.floor(y);
		const j1 = j0 + 1;
		const k0 = Math.floor(z);
		const k1 = k0 + 1;
		const fx = x - i0;
		const fy = y - j0;
		const fz = z - k0;
		const u = smoothStep(fx);
		const v = smoothStep(fy);
		const w = smoothStep(fz);
		const p = this.p;
		const A0 = p[i0 & 255] + (j0 & 255),
			AA0 = p[A0] + (k0 & 255),
			AB0 = p[A0 + 1] + (k0 & 255);
		const B0 = p[i1 & 255] + (j0 & 255),
			BA0 = p[B0] + (k0 & 255),
			BB0 = p[B0 + 1] + (k0 & 255);
		const A1 = p[i0 & 255] + (j1 & 255),
			AA1 = p[A1] + (k1 & 255),
			AB1 = p[A1 + 1] + (k1 & 255);
		const B1 = p[i1 & 255] + (j1 & 255),
			BA1 = p[B1] + (k1 & 255),
			BB1 = p[B1 + 1] + (k1 & 255);
		const v000 = p[AA0] / 255;
		const v100 = p[BA0] / 255;
		const v010 = p[AB0] / 255;
		const v110 = p[BB0] / 255;
		const v001 = p[AA1] / 255;
		const v101 = p[BA1] / 255;
		const v011 = p[AB1] / 255;
		const v111 = p[BB1] / 255;
		const mixXY0 = interpolate(v, interpolate(u, v000, v100), interpolate(u, v010, v110));
		const mixXY1 = interpolate(v, interpolate(u, v001, v101), interpolate(u, v011, v111));
		return interpolate(w, mixXY0, mixXY1);
	};

	ForgeNoise.prototype.getCellPoint2D = function(i, j) {
		const p = this.p;
		const ii = i & 255;
		const jj = j & 255;
		const A = p[jj];
		const hash1 = p[(ii + A) & 255];
		const hash2 = p[(ii + 1 + A) & 255];
		return [i + hash1 / 256, j + hash2 / 256];
	};

	ForgeNoise.prototype.getCellPoint3D = function(i, j, k) {
		const p = this.p;
		const ii = i & 255;
		const jj = j & 255;
		const kk = k & 255;
		const A = p[jj + p[kk]];
		const hash1 = p[(ii + A) & 255];
		const hash2 = p[(ii + 1 + A) & 255];
		const hash3 = p[(ii + 2 + A) & 255];
		return [i + hash1 / 256, j + hash2 / 256, k + hash3 / 256];
	};

	ForgeNoise.prototype.generateWorley2D = function(x, y) {
		const i0 = Math.floor(x);
		const j0 = Math.floor(y);
		let minDist = Infinity;
		for (let di = -1; di <= 1; di++) {
			for (let dj = -1; dj <= 1; dj++) {
				const i = i0 + di;
				const j = j0 + dj;
				const [px, py] = this.getCellPoint2D(i, j);
				const dx = px - x;
				const dy = py - y;
				const dist = Math.sqrt(dx * dx + dy * dy);
				minDist = Math.min(minDist, dist);
			}
		}
		return minDist;
	};

	ForgeNoise.prototype.generateWorley3D = function(x, y, z) {
		const i0 = Math.floor(x);
		const j0 = Math.floor(y);
		const k0 = Math.floor(z);
		let minDist = Infinity;
		for (let di = -1; di <= 1; di++) {
			for (let dj = -1; dj <= 1; dj++) {
				for (let dk = -1; dk <= 1; dk++) {
					const i = i0 + di;
					const j = j0 + dj;
					const k = k0 + dk;
					const [px, py, pz] = this.getCellPoint3D(i, j, k);
					const dx = px - x;
					const dy = py - y;
					const dz = pz - z;
					const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
					minDist = Math.min(minDist, dist);
				}
			}
		}
		return minDist;
	};

	ForgeNoise.prototype.generateVoronoi2D = function(x, y) {
		const i0 = Math.floor(x);
		const j0 = Math.floor(y);
		let minDist = Infinity;
		let closestI, closestJ;
		for (let di = -1; di <= 1; di++) {
			for (let dj = -1; dj <= 1; dj++) {
				const i = i0 + di;
				const j = j0 + dj;
				const [px, py] = this.getCellPoint2D(i, j);
				const dx = px - x;
				const dy = py - y;
				const dist = dx * dx + dy * dy;
				if (dist < minDist) {
					minDist = dist;
					closestI = i;
					closestJ = j;
				}
			}
		}
		return {
			distance: Math.sqrt(minDist),
			cellI: closestI,
			cellJ: closestJ
		};
	};

	ForgeNoise.prototype.generateFractal2D = function(x, y, options = {}) {
		const octaves = options.octaves || 4;
		const persistence = options.persistence || 0.5;
		const lacunarity = options.lacunarity || 2.0;
		const scale = options.scale || 1.0;
		const normalized = options.normalized !== false;
		const ridged = options.ridged || false;
		const turbulence = options.turbulence || false;

		let total = 0;
		let frequency = scale;
		let amplitude = 1;
		let maxValue = 0;

		for (let i = 0; i < octaves; i++) {
			let noiseValue = this.generate2D(x * frequency, y * frequency);
			if (ridged) noiseValue = 1 - Math.abs(noiseValue);
			else if (turbulence) noiseValue = Math.abs(noiseValue);
			total += noiseValue * amplitude;
			maxValue += amplitude * (ridged ? 1 : 1);
			amplitude *= persistence;
			frequency *= lacunarity;
		}

		return normalized ? total / maxValue : total;
	};

	ForgeNoise.prototype.generateFractal3D = function(x, y, z, options = {}) {
		const octaves = options.octaves || 4;
		const persistence = options.persistence || 0.5;
		const lacunarity = options.lacunarity || 2.0;
		const scale = options.scale || 1.0;
		const normalized = options.normalized !== false;
		const ridged = options.ridged || false;
		const turbulence = options.turbulence || false;

		let total = 0;
		let frequency = scale;
		let amplitude = 1;
		let maxValue = 0;

		for (let i = 0; i < octaves; i++) {
			let noiseValue = this.generate3D(x * frequency, y * frequency, z * frequency);
			if (ridged) noiseValue = 1 - Math.abs(noiseValue);
			else if (turbulence) noiseValue = Math.abs(noiseValue);
			total += noiseValue * amplitude;
			maxValue += amplitude * (ridged ? 1 : 1);
			amplitude *= persistence;
			frequency *= lacunarity;
		}

		return normalized ? total / maxValue : total;
	};

	ForgeNoise.prototype.generateSimplexFractal2D = function(x, y, options = {}) {
		const octaves = options.octaves || 4;
		const persistence = options.persistence || 0.5;
		const lacunarity = options.lacunarity || 2.0;
		const scale = options.scale || 1.0;
		const normalized = options.normalized !== false;
		const ridged = options.ridged || false;
		const turbulence = options.turbulence || false;

		let total = 0;
		let frequency = scale;
		let amplitude = 1;
		let maxValue = 0;

		for (let i = 0; i < octaves; i++) {
			let noiseValue = this.generateSimplex2D(x * frequency, y * frequency);
			if (ridged) noiseValue = 1 - Math.abs(noiseValue);
			else if (turbulence) noiseValue = Math.abs(noiseValue);
			total += noiseValue * amplitude;
			maxValue += amplitude * (ridged ? 1 : 1);
			amplitude *= persistence;
			frequency *= lacunarity;
		}

		return normalized ? total / maxValue : total;
	};

	ForgeNoise.prototype.generateSimplexFractal3D = function(x, y, z, options = {}) {
		const octaves = options.octaves || 4;
		const persistence = options.persistence || 0.5;
		const lacunarity = options.lacunarity || 2.0;
		const scale = options.scale || 1.0;
		const normalized = options.normalized !== false;
		const ridged = options.ridged || false;
		const turbulence = options.turbulence || false;

		let total = 0;
		let frequency = scale;
		let amplitude = 1;
		let maxValue = 0;

		for (let i = 0; i < octaves; i++) {
			let noiseValue = this.generateSimplex3D(x * frequency, y * frequency, z * frequency);
			if (ridged) noiseValue = 1 - Math.abs(noiseValue);
			else if (turbulence) noiseValue = Math.abs(noiseValue);
			total += noiseValue * amplitude;
			maxValue += amplitude * (ridged ? 1 : 1);
			amplitude *= persistence;
			frequency *= lacunarity;
		}

		return normalized ? total / maxValue : total;
	};

	ForgeNoise.prototype.generate4D = function(x, y, z, w) {
		const X = Math.floor(x) & 255;
		const Y = Math.floor(y) & 255;
		const Z = Math.floor(z) & 255;
		const W = Math.floor(w) & 255;
		x -= Math.floor(x);
		y -= Math.floor(y);
		z -= Math.floor(z);
		w -= Math.floor(w);
		const u = smoothStep(x);
		const v = smoothStep(y);
		const w1 = smoothStep(z);
		const t = smoothStep(w);
		const p = this.p;

		const A = p[X] + Y,
			AA = p[A] + Z,
			AAA = p[AA] + W,
			AAB = p[AA + 1] + W;
		const B = p[X + 1] + Y,
			BA = p[B] + Z,
			BAA = p[BA] + W,
			BAB = p[BA + 1] + W;
		const AB = p[A + 1] + Z,
			ABA = p[AB] + W,
			ABB = p[AB + 1] + W;
		const BB = p[B + 1] + Z,
			BBA = p[BB] + W,
			BBB = p[BB + 1] + W;

		const hash = (h) => h & 31;
		const grad4D = (hash, x, y, z, w) => {
			const g = grad4[hash];
			return dot4D(g, x, y, z, w);
		};

		return interpolate(t,
			interpolate(w1,
				interpolate(v,
					interpolate(u, grad4D(p[AAA], x, y, z, w), grad4D(p[BAA], x - 1, y, z, w)),
					interpolate(u, grad4D(p[ABA], x, y - 1, z, w), grad4D(p[BBA], x - 1, y - 1, z, w))
				),
				interpolate(v,
					interpolate(u, grad4D(p[AAB], x, y, z - 1, w), grad4D(p[BAB], x - 1, y, z - 1, w)),
					interpolate(u, grad4D(p[ABB], x, y - 1, z - 1, w), grad4D(p[BBB], x - 1, y - 1, z - 1, w))
				)
			),
			interpolate(w1,
				interpolate(v,
					interpolate(u, grad4D(p[AAA + 1], x, y, z, w - 1), grad4D(p[BAA + 1], x - 1, y, z, w - 1)),
					interpolate(u, grad4D(p[ABA + 1], x, y - 1, z, w - 1), grad4D(p[BBA + 1], x - 1, y - 1, z, w - 1))
				),
				interpolate(v,
					interpolate(u, grad4D(p[AAB + 1], x, y, z - 1, w - 1), grad4D(p[BAB + 1], x - 1, y, z - 1, w - 1)),
					interpolate(u, grad4D(p[ABB + 1], x, y - 1, z - 1, w - 1), grad4D(p[BBB + 1], x - 1, y - 1, z - 1, w - 1))
				)
			)
		);
	};

	ForgeNoise.prototype.generateSimplex4D = function(x, y, z, w) {
		const F4 = (Math.sqrt(5) - 1) / 4;
		const G4 = (5 - Math.sqrt(5)) / 20;
		const s = (x + y + z + w) * F4;
		const i = Math.floor(x + s);
		const j = Math.floor(y + s);
		const k = Math.floor(z + s);
		const l = Math.floor(w + s);
		const t = (i + j + k + l) * G4;
		const X0 = i - t;
		const Y0 = j - t;
		const Z0 = k - t;
		const W0 = l - t;
		const x0 = x - X0;
		const y0 = y - Y0;
		const z0 = z - Z0;
		const w0 = w - W0;

		let c1 = (x0 > y0) ? 32 : 0;
		let c2 = (x0 > z0) ? 16 : 0;
		let c3 = (y0 > z0) ? 8 : 0;
		let c4 = (x0 > w0) ? 4 : 0;
		let c5 = (y0 > w0) ? 2 : 0;
		let c6 = (z0 > w0) ? 1 : 0;
		const c = c1 + c2 + c3 + c4 + c5 + c6;
		const [i1, j1, k1, l1, i2, j2, k2, l2, i3, j3, k3, l3] = [
			[1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0],
			[1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0],
			[0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0],
			[1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0],
			[0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
			[0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0],
			[0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0],
			[0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0],
			[0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0],
			[0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0],
			[0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0],
			[1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1],
			[0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1],
			[0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1],
			[0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1],
			[0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1]
		][c >>> 0];

		const x1 = x0 - i1 + G4,
			y1 = y0 - j1 + G4,
			z1 = z0 - k1 + G4,
			w1 = w0 - l1 + G4;
		const x2 = x0 - i2 + 2 * G4,
			y2 = y0 - j2 + 2 * G4,
			z2 = z0 - k2 + 2 * G4,
			w2 = w0 - l2 + 2 * G4;
		const x3 = x0 - i3 + 3 * G4,
			y3 = y0 - j3 + 3 * G4,
			z3 = z0 - k3 + 3 * G4,
			w3 = w0 - l3 + 3 * G4;
		const x4 = x0 - 1 + 4 * G4,
			y4 = y0 - 1 + 4 * G4,
			z4 = z0 - 1 + 4 * G4,
			w4 = w0 - 1 + 4 * G4;

		const ii = i & 255,
			jj = j & 255,
			kk = k & 255,
			ll = l & 255;
		const p = this.p;
		const gi0 = p[ii + p[jj + p[kk + p[ll]]]] % 32;
		const gi1 = p[ii + i1 + p[jj + j1 + p[kk + k1 + p[ll + l1]]]] % 32;
		const gi2 = p[ii + i2 + p[jj + j2 + p[kk + k2 + p[ll + l2]]]] % 32;
		const gi3 = p[ii + i3 + p[jj + j3 + p[kk + k3 + p[ll + l3]]]] % 32;
		const gi4 = p[ii + 1 + p[jj + 1 + p[kk + 1 + p[ll + 1]]]] % 32;

		let n0 = 0,
			n1 = 0,
			n2 = 0,
			n3 = 0,
			n4 = 0;
		let t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0 - w0 * w0;
		if (t0 >= 0) {
			t0 *= t0;
			n0 = t0 * t0 * dot4D(grad4[gi0], x0, y0, z0, w0);
		}
		let t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1 - w1 * w1;
		if (t1 >= 0) {
			t1 *= t1;
			n1 = t1 * t1 * dot4D(grad4[gi1], x1, y1, z1, w1);
		}
		let t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2 - w2 * w2;
		if (t2 >= 0) {
			t2 *= t2;
			n2 = t2 * t2 * dot4D(grad4[gi2], x2, y2, z2, w2);
		}
		let t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3 - w3 * w3;
		if (t3 >= 0) {
			t3 *= t3;
			n3 = t3 * t3 * dot4D(grad4[gi3], x3, y3, z3, w3);
		}
		let t4 = 0.6 - x4 * x4 - y4 * y4 - z4 * z4 - w4 * w4;
		if (t4 >= 0) {
			t4 *= t4;
			n4 = t4 * t4 * dot4D(grad4[gi4], x4, y4, z4, w4);
		}
		return 27.0 * (n0 + n1 + n2 + n3 + n4);
	};

	ForgeNoise.prototype.generateTiling2D = function(x, y, options = {}) {
		const periodX = options.periodX || 256;
		const periodY = options.periodY || 256;
		const xWrapped = x % periodX;
		const yWrapped = y % periodY;
		return this.generate2D(xWrapped, yWrapped);
	};

	ForgeNoise.prototype.warp2D = function(x, y, options = {}) {
		const warpStrength = options.warpStrength || 1.0;
		const warpScale = options.warpScale || 1.0;
		const iterations = options.iterations || 1;
		let qx = x;
		let qy = y;

		for (let i = 0; i < iterations; i++) {
			const offsetX = this.generate2D(qx * warpScale, qy * warpScale);
			const offsetY = this.generate2D(qx * warpScale + 5.2, qy * warpScale + 1.3);

			qx += offsetX * warpStrength;
			qy += offsetY * warpStrength;
		}

		const finalNoise = this.generate2D(qx, qy);
		return finalNoise;
	};

	ForgeNoise.prototype.generate2D01 = function(x, y) {
		return (this.generate2D(x, y) + 1) * 0.5;
	};
	ForgeNoise.prototype.generate3D01 = function(x, y, z) {
		return (this.generate3D(x, y, z) + 1) * 0.5;
	};
	ForgeNoise.prototype.generateSimplex2D01 = function(x, y) {
		return (this.generateSimplex2D(x, y) + 1) * 0.5;
	};
	ForgeNoise.prototype.generateSimplex3D01 = function(x, y, z) {
		return (this.generateSimplex3D(x, y, z) + 1) * 0.5;
	};
	ForgeNoise.prototype.generateValue2D01 = function(x, y) {
		return this.generateValue2D(x, y);
	};
	ForgeNoise.prototype.generateValue3D01 = function(x, y, z) {
		return this.generateValue3D(x, y, z);
	};
	ForgeNoise.prototype.generateWorley2D01 = function(x, y) {
		return this.generateWorley2D(x, y) / 1.414;
	};
	ForgeNoise.prototype.generateWorley3D01 = function(x, y, z) {
		return this.generateWorley3D(x, y, z) / 1.732;
	};
	ForgeNoise.prototype.generate4D01 = function(x, y, z, w) {
		return (this.generate4D(x, y, z, w) + 1) * 0.5;
	};
	ForgeNoise.prototype.generateSimplex4D01 = function(x, y, z, w) {
		return (this.generateSimplex4D(x, y, z, w) + 1) * 0.5;
	};
	ForgeNoise.prototype.generateTiling2D01 = function(x, y, options) {
		return (this.generateTiling2D(x, y, options) + 1) * 0.5;
	};

	if (typeof module !== 'undefined' && module.exports) {
		module.exports = ForgeNoise;
	} else {
		root.ForgeNoise = ForgeNoise;
	}
})(this);
