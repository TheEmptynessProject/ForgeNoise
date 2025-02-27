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

	const smoothStep = t => t * t * t * (t * (t * 6 - 15) + 10);
	const lerp = (t, a, b) => a + t * (b - a);
	const dot2D = (g, x, y) => g[0] * x + g[1] * y;
	const dot3D = (g, x, y, z) => g[0] * x + g[1] * y + g[2] * z;
	const dot4D = (g, x, y, z, w) => g[0] * x + g[1] * y + g[2] * z + g[3] * w;

	function gradient2D(hash, x, y) {
		let h = hash & 7;

		let useY = (h >> 2); 
		let u = useY ? y : x;
		let v = useY ? x : y;
		return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
	}

	function gradient3D(hash, x, y, z) {
		let h = hash & 15;
		let u = (h < 8) ? x : y;
		let v = (h < 4) ? y : ((h === 12 || h === 14) ? x : z);
		return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
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

	ForgeNoise.prototype.generate2D = function(x, y) {
		let ix = Math.floor(x), iy = Math.floor(y);
		let X = ix & 255, Y = iy & 255;
		let fx = x - ix, fy = y - iy;
		let u = smoothStep(fx);
		let v = smoothStep(fy);
		let p = this.p;
		let A = p[X] + Y, B = p[X + 1] + Y;
		let gA = gradient2D(p[A], fx, fy);
		let gB = gradient2D(p[B], fx - 1, fy);
		let gA1 = gradient2D(p[A + 1], fx, fy - 1);
		let gB1 = gradient2D(p[B + 1], fx - 1, fy - 1);
		return lerp(v, lerp(u, gA, gB), lerp(u, gA1, gB1));
	};

	ForgeNoise.prototype.generate3D = function(x, y, z) {
		let ix = Math.floor(x), iy = Math.floor(y), iz = Math.floor(z);
		let X = ix & 255, Y = iy & 255, Z = iz & 255;
		let fx = x - ix, fy = y - iy, fz = z - iz;
		let u = smoothStep(fx), v = smoothStep(fy), w = smoothStep(fz);
		let p = this.p;
		let A = p[X] + Y, B = p[X + 1] + Y;
		let AA = p[A] + Z, AB = p[A + 1] + Z;
		let BA = p[B] + Z, BB = p[B + 1] + Z;
		return lerp(w,
			lerp(v,
				lerp(u, gradient3D(p[AA], fx, fy, fz), gradient3D(p[BA], fx - 1, fy, fz)),
				lerp(u, gradient3D(p[AB], fx, fy - 1, fz), gradient3D(p[BB], fx - 1, fy - 1, fz))
			),
			lerp(v,
				lerp(u, gradient3D(p[AA + 1], fx, fy, fz - 1), gradient3D(p[BA + 1], fx - 1, fy, fz - 1)),
				lerp(u, gradient3D(p[AB + 1], fx, fy - 1, fz - 1), gradient3D(p[BB + 1], fx - 1, fy - 1, fz - 1))
			)
		);
	};

	ForgeNoise.prototype.generateSimplex2D = function(x, y) {
		const F2 = 0.5 * (Math.sqrt(3) - 1);
		const G2 = (3 - Math.sqrt(3)) / 6;
		let s = (x + y) * F2;
		let i = Math.floor(x + s), j = Math.floor(y + s);
		let t = (i + j) * G2;
		let X0 = i - t, Y0 = j - t;
		let x0 = x - X0, y0 = y - Y0;
		let i1, j1;
		if (x0 > y0) { i1 = 1; j1 = 0; } else { i1 = 0; j1 = 1; }
		let x1 = x0 - i1 + G2, y1 = y0 - j1 + G2;
		let x2 = x0 - 1 + 2 * G2, y2 = y0 - 1 + 2 * G2;
		let ii = i & 255, jj = j & 255;
		let p = this.p;
		let gi0 = p[ii + p[jj]] % 12;
		let gi1 = p[ii + i1 + p[jj + j1]] % 12;
		let gi2 = p[ii + 1 + p[jj + 1]] % 12;
		let n0 = 0, n1 = 0, n2 = 0;
		let t0 = 0.5 - x0 * x0 - y0 * y0;
		if (t0 >= 0) { t0 *= t0; n0 = t0 * t0 * dot2D(grad2[gi0], x0, y0); }
		let t1 = 0.5 - x1 * x1 - y1 * y1;
		if (t1 >= 0) { t1 *= t1; n1 = t1 * t1 * dot2D(grad2[gi1], x1, y1); }
		let t2 = 0.5 - x2 * x2 - y2 * y2;
		if (t2 >= 0) { t2 *= t2; n2 = t2 * t2 * dot2D(grad2[gi2], x2, y2); }
		return 70 * (n0 + n1 + n2);
	};

	ForgeNoise.prototype.generateSimplex3D = function(x, y, z) {
		const F3 = 1 / 3, G3 = 1 / 6;
		let s = (x + y + z) * F3;
		let i = Math.floor(x + s), j = Math.floor(y + s), k = Math.floor(z + s);
		let t = (i + j + k) * G3;
		let X0 = i - t, Y0 = j - t, Z0 = k - t;
		let x0 = x - X0, y0 = y - Y0, z0 = z - Z0;
		let i1, j1, k1, i2, j2, k2;
		if (x0 >= y0) {
			if (y0 >= z0) { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0; }
			else if (x0 >= z0) { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1; }
			else { i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1; }
		} else {
			if (y0 < z0) { i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1; }
			else if (x0 < z0) { i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1; }
			else { i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0; }
		}
		let x1 = x0 - i1 + G3, y1 = y0 - j1 + G3, z1 = z0 - k1 + G3;
		let x2 = x0 - i2 + 2 * G3, y2 = y0 - j2 + 2 * G3, z2 = z0 - k2 + 2 * G3;
		let x3 = x0 - 1 + 3 * G3, y3 = y0 - 1 + 3 * G3, z3 = z0 - 1 + 3 * G3;
		let ii = i & 255, jj = j & 255, kk = k & 255;
		let p = this.p;
		let gi0 = p[ii + p[jj + p[kk]]] % 12;
		let gi1 = p[ii + i1 + p[jj + j1 + p[kk + k1]]] % 12;
		let gi2 = p[ii + i2 + p[jj + j2 + p[kk + k2]]] % 12;
		let gi3 = p[ii + 1 + p[jj + 1 + p[kk + 1]]] % 12;
		let n0 = 0, n1 = 0, n2 = 0, n3 = 0;
		let t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
		if (t0 >= 0) { t0 *= t0; n0 = t0 * t0 * dot3D(grad3[gi0], x0, y0, z0); }
		let t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
		if (t1 >= 0) { t1 *= t1; n1 = t1 * t1 * dot3D(grad3[gi1], x1, y1, z1); }
		let t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
		if (t2 >= 0) { t2 *= t2; n2 = t2 * t2 * dot3D(grad3[gi2], x2, y2, z2); }
		let t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
		if (t3 >= 0) { t3 *= t3; n3 = t3 * t3 * dot3D(grad3[gi3], x3, y3, z3); }
		return 32 * (n0 + n1 + n2 + n3);
	};

	ForgeNoise.prototype.generateValue2D = function(x, y) {
		let i0 = Math.floor(x), i1 = i0 + 1;
		let j0 = Math.floor(y), j1 = j0 + 1;
		let fx = x - i0, fy = y - j0;
		let u = smoothStep(fx), v = smoothStep(fy);
		let p = this.p;
		let a00 = p[(i0 + p[j0 & 255]) & 255] / 255;
		let a01 = p[(i0 + p[j1 & 255]) & 255] / 255;
		let a10 = p[(i1 + p[j0 & 255]) & 255] / 255;
		let a11 = p[(i1 + p[j1 & 255]) & 255] / 255;
		let mixX0 = lerp(u, a00, a10);
		let mixX1 = lerp(u, a01, a11);
		return lerp(v, mixX0, mixX1);
	};

	ForgeNoise.prototype.generateValue3D = function(x, y, z) {
		let i0 = Math.floor(x), i1 = i0 + 1;
		let j0 = Math.floor(y), j1 = j0 + 1;
		let k0 = Math.floor(z), k1 = k0 + 1;
		let fx = x - i0, fy = y - j0, fz = z - k0;
		let u = smoothStep(fx), v = smoothStep(fy), w = smoothStep(fz);
		let p = this.p;
		let A0 = p[i0 & 255] + (j0 & 255),
			AA0 = p[A0] + (k0 & 255),
			AB0 = p[A0 + 1] + (k0 & 255);
		let B0 = p[i1 & 255] + (j0 & 255),
			BA0 = p[B0] + (k0 & 255),
			BB0 = p[B0 + 1] + (k0 & 255);
		let A1 = p[i0 & 255] + (j1 & 255),
			AA1 = p[A1] + (k1 & 255),
			AB1 = p[A1 + 1] + (k1 & 255);
		let B1 = p[i1 & 255] + (j1 & 255),
			BA1 = p[B1] + (k1 & 255),
			BB1 = p[B1 + 1] + (k1 & 255);
		let v000 = p[AA0] / 255;
		let v100 = p[BA0] / 255;
		let v010 = p[AB0] / 255;
		let v110 = p[BB0] / 255;
		let v001 = p[AA1] / 255;
		let v101 = p[BA1] / 255;
		let v011 = p[AB1] / 255;
		let v111 = p[BB1] / 255;
		let mixXY0 = lerp(v, lerp(u, v000, v100), lerp(u, v010, v110));
		let mixXY1 = lerp(v, lerp(u, v001, v101), lerp(u, v011, v111));
		return lerp(w, mixXY0, mixXY1);
	};

	ForgeNoise.prototype.getCellPoint2D = function(i, j) {
		let p = this.p;
		let ii = i & 255, jj = j & 255;
		let A = p[jj];
		let hash1 = p[(ii + A) & 255];
		let hash2 = p[(ii + 1 + A) & 255];
		return [i + hash1 / 256, j + hash2 / 256];
	};

	ForgeNoise.prototype.getCellPoint3D = function(i, j, k) {
		let p = this.p;
		let ii = i & 255, jj = j & 255, kk = k & 255;
		let A = p[jj + p[kk]];
		let hash1 = p[(ii + A) & 255];
		let hash2 = p[(ii + 1 + A) & 255];
		let hash3 = p[(ii + 2 + A) & 255];
		return [i + hash1 / 256, j + hash2 / 256, k + hash3 / 256];
	};

	ForgeNoise.prototype.generateWorley2D = function(x, y) {
		let i0 = Math.floor(x), j0 = Math.floor(y);
		let minDist = Infinity;
		const offsets = [-1, 0, 1];
		for (let di of offsets) {
			for (let dj of offsets) {
				let cell = this.getCellPoint2D(i0 + di, j0 + dj);
				let dx = cell[0] - x, dy = cell[1] - y;
				let d = Math.hypot(dx, dy);
				if (d < minDist) minDist = d;
			}
		}
		return minDist;
	};

	ForgeNoise.prototype.generateWorley3D = function(x, y, z) {
		let i0 = Math.floor(x), j0 = Math.floor(y), k0 = Math.floor(z);
		let minDist = Infinity;
		const offsets = [-1, 0, 1];
		for (let di of offsets) {
			for (let dj of offsets) {
				for (let dk of offsets) {
					let cell = this.getCellPoint3D(i0 + di, j0 + dj, k0 + dk);
					let dx = cell[0] - x, dy = cell[1] - y, dz = cell[2] - z;
					let d = Math.sqrt(dx * dx + dy * dy + dz * dz);
					if (d < minDist) minDist = d;
				}
			}
		}
		return minDist;
	};

	ForgeNoise.prototype.generateVoronoi2D = function(x, y) {
		let i0 = Math.floor(x), j0 = Math.floor(y);
		let minDist = Infinity, closestI, closestJ;
		const offsets = [-1, 0, 1];
		for (let di of offsets) {
			for (let dj of offsets) {
				let i = i0 + di, j = j0 + dj;
				let cell = this.getCellPoint2D(i, j);
				let dx = cell[0] - x, dy = cell[1] - y;
				let d2 = dx * dx + dy * dy;
				if (d2 < minDist) { minDist = d2; closestI = i; closestJ = j; }
			}
		}
		return { distance: Math.sqrt(minDist), cellI: closestI, cellJ: closestJ };
	};

	ForgeNoise.prototype.generateFractal2D = function(x, y, options = {}) {
		let octaves = options.octaves || 4;
		let persistence = options.persistence || 0.5;
		let lacunarity = options.lacunarity || 2.0;
		let scale = options.scale || 1.0;
		let normalized = options.normalized !== false;
		let ridged = options.ridged || false;
		let turbulence = options.turbulence || false;
		let total = 0, frequency = scale, amplitude = 1, maxValue = 0;
		for (let i = 0; i < octaves; i++) {
			let noiseValue = this.generate2D(x * frequency, y * frequency);
			if (ridged) noiseValue = 1 - Math.abs(noiseValue);
			else if (turbulence) noiseValue = Math.abs(noiseValue);
			total += noiseValue * amplitude;
			maxValue += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		return normalized ? total / maxValue : total;
	};

	ForgeNoise.prototype.generateFractal3D = function(x, y, z, options = {}) {
		let octaves = options.octaves || 4;
		let persistence = options.persistence || 0.5;
		let lacunarity = options.lacunarity || 2.0;
		let scale = options.scale || 1.0;
		let normalized = options.normalized !== false;
		let ridged = options.ridged || false;
		let turbulence = options.turbulence || false;
		let total = 0, frequency = scale, amplitude = 1, maxValue = 0;
		for (let i = 0; i < octaves; i++) {
			let noiseValue = this.generate3D(x * frequency, y * frequency, z * frequency);
			if (ridged) noiseValue = 1 - Math.abs(noiseValue);
			else if (turbulence) noiseValue = Math.abs(noiseValue);
			total += noiseValue * amplitude;
			maxValue += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		return normalized ? total / maxValue : total;
	};

	ForgeNoise.prototype.generateSimplexFractal2D = function(x, y, options = {}) {
		let octaves = options.octaves || 4;
		let persistence = options.persistence || 0.5;
		let lacunarity = options.lacunarity || 2.0;
		let scale = options.scale || 1.0;
		let normalized = options.normalized !== false;
		let ridged = options.ridged || false;
		let turbulence = options.turbulence || false;
		let total = 0, frequency = scale, amplitude = 1, maxValue = 0;
		for (let i = 0; i < octaves; i++) {
			let noiseValue = this.generateSimplex2D(x * frequency, y * frequency);
			if (ridged) noiseValue = 1 - Math.abs(noiseValue);
			else if (turbulence) noiseValue = Math.abs(noiseValue);
			total += noiseValue * amplitude;
			maxValue += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		return normalized ? total / maxValue : total;
	};

	ForgeNoise.prototype.generateSimplexFractal3D = function(x, y, z, options = {}) {
		let octaves = options.octaves || 4;
		let persistence = options.persistence || 0.5;
		let lacunarity = options.lacunarity || 2.0;
		let scale = options.scale || 1.0;
		let normalized = options.normalized !== false;
		let ridged = options.ridged || false;
		let turbulence = options.turbulence || false;
		let total = 0, frequency = scale, amplitude = 1, maxValue = 0;
		for (let i = 0; i < octaves; i++) {
			let noiseValue = this.generateSimplex3D(x * frequency, y * frequency, z * frequency);
			if (ridged) noiseValue = 1 - Math.abs(noiseValue);
			else if (turbulence) noiseValue = Math.abs(noiseValue);
			total += noiseValue * amplitude;
			maxValue += amplitude;
			amplitude *= persistence;
			frequency *= lacunarity;
		}
		return normalized ? total / maxValue : total;
	};

	ForgeNoise.prototype.generate4D = function(x, y, z, w) {
		let ix = Math.floor(x), iy = Math.floor(y), iz = Math.floor(z), iw = Math.floor(w);
		let X = ix & 255, Y = iy & 255, Z = iz & 255, W = iw & 255;
		let fx = x - ix, fy = y - iy, fz = z - iz, fw = w - iw;
		let u = smoothStep(fx), v = smoothStep(fy), w1 = smoothStep(fz), t = smoothStep(fw);
		let p = this.p;
		let A = p[X] + Y, AA = p[A] + Z, AAA = p[AA] + W, AAB = p[AA + 1] + W;
		let B = p[X + 1] + Y, BA = p[B] + Z, BAA = p[BA] + W, BAB = p[BA + 1] + W;
		let AB = p[A + 1] + Z, ABA = p[AB] + W, ABB = p[AB + 1] + W;
		let BB = p[B + 1] + Z, BBA = p[BB] + W, BBB = p[BB + 1] + W;

		let grad4D = (hash, x, y, z, w) => dot4D(grad4[hash & 31], x, y, z, w);
		return lerp(t,
			lerp(w1,
				lerp(v,
					lerp(u, grad4D(p[AAA], fx, fy, fz, fw), grad4D(p[BAA], fx - 1, fy, fz, fw)),
					lerp(u, grad4D(p[ABA], fx, fy - 1, fz, fw), grad4D(p[BBA], fx - 1, fy - 1, fz, fw))
				),
				lerp(v,
					lerp(u, grad4D(p[AAB], fx, fy, fz - 1, fw), grad4D(p[BAB], fx - 1, fy, fz - 1, fw)),
					lerp(u, grad4D(p[ABB], fx, fy - 1, fz - 1, fw), grad4D(p[BBB], fx - 1, fy - 1, fz - 1, fw))
				)
			),
			lerp(w1,
				lerp(v,
					lerp(u, grad4D(p[AAA + 1], fx, fy, fz, fw - 1), grad4D(p[BAA + 1], fx - 1, fy, fz, fw - 1)),
					lerp(u, grad4D(p[ABA + 1], fx, fy - 1, fz, fw - 1), grad4D(p[BBA + 1], fx - 1, fy - 1, fz, fw - 1))
				),
				lerp(v,
					lerp(u, grad4D(p[AAB + 1], fx, fy, fz - 1, fw - 1), grad4D(p[BAB + 1], fx - 1, fy, fz - 1, fw - 1)),
					lerp(u, grad4D(p[ABB + 1], fx, fy - 1, fz - 1, fw - 1), grad4D(p[BBB + 1], fx - 1, fy - 1, fz - 1, fw - 1))
				)
			)
		);
	};

	ForgeNoise.prototype.generateSimplex4D = function(x, y, z, w) {
		const F4 = (Math.sqrt(5) - 1) / 4;
		const G4 = (5 - Math.sqrt(5)) / 20;
		let s = (x + y + z + w) * F4;
		let i = Math.floor(x + s), j = Math.floor(y + s), k = Math.floor(z + s), l = Math.floor(w + s);
		let t = (i + j + k + l) * G4;
		let X0 = i - t, Y0 = j - t, Z0 = k - t, W0 = l - t;
		let x0 = x - X0, y0 = y - Y0, z0 = z - Z0, w0 = w - W0;
		let c1 = (x0 > y0) ? 32 : 0;
		let c2 = (x0 > z0) ? 16 : 0;
		let c3 = (y0 > z0) ? 8 : 0;
		let c4 = (x0 > w0) ? 4 : 0;
		let c5 = (y0 > w0) ? 2 : 0;
		let c6 = (z0 > w0) ? 1 : 0;
		const c = c1 + c2 + c3 + c4 + c5 + c6;

		const lookup = [
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
		][c];
		let i1 = lookup[0], j1 = lookup[1], k1 = lookup[2], l1 = lookup[3];
		let i2 = lookup[4], j2 = lookup[5], k2 = lookup[6], l2 = lookup[7];
		let i3 = lookup[8], j3 = lookup[9], k3 = lookup[10], l3 = lookup[11];
		let x1 = x0 - i1 + G4, y1 = y0 - j1 + G4, z1 = z0 - k1 + G4, w1 = w0 - l1 + G4;
		let x2 = x0 - i2 + 2 * G4, y2 = y0 - j2 + 2 * G4, z2 = z0 - k2 + 2 * G4, w2 = w0 - l2 + 2 * G4;
		let x3 = x0 - i3 + 3 * G4, y3 = y0 - j3 + 3 * G4, z3 = z0 - k3 + 3 * G4, w3 = w0 - l3 + 3 * G4;
		let x4 = x0 - 1 + 4 * G4, y4 = y0 - 1 + 4 * G4, z4 = z0 - 1 + 4 * G4, w4 = w0 - 1 + 4 * G4;
		let ii = i & 255, jj = j & 255, kk = k & 255, ll = l & 255;
		let p = this.p;
		let gi0 = p[ii + p[jj + p[kk + p[ll]]]] % 32;
		let gi1 = p[ii + i1 + p[jj + j1 + p[kk + k1 + p[ll + l1]]]] % 32;
		let gi2 = p[ii + i2 + p[jj + j2 + p[kk + k2 + p[ll + l2]]]] % 32;
		let gi3 = p[ii + i3 + p[jj + j3 + p[kk + k3 + p[ll + l3]]]] % 32;
		let gi4 = p[ii + 1 + p[jj + 1 + p[kk + 1 + p[ll + 1]]]] % 32;
		let n0 = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0;
		let t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0 - w0 * w0;
		if (t0 >= 0) { t0 *= t0; n0 = t0 * t0 * dot4D(grad4[gi0], x0, y0, z0, w0); }
		let t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1 - w1 * w1;
		if (t1 >= 0) { t1 *= t1; n1 = t1 * t1 * dot4D(grad4[gi1], x1, y1, z1, w1); }
		let t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2 - w2 * w2;
		if (t2 >= 0) { t2 *= t2; n2 = t2 * t2 * dot4D(grad4[gi2], x2, y2, z2, w2); }
		let t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3 - w3 * w3;
		if (t3 >= 0) { t3 *= t3; n3 = t3 * t3 * dot4D(grad4[gi3], x3, y3, z3, w3); }
		let t4 = 0.6 - x4 * x4 - y4 * y4 - z4 * z4 - w4 * w4;
		if (t4 >= 0) { t4 *= t4; n4 = t4 * t4 * dot4D(grad4[gi4], x4, y4, z4, w4); }
		return 27 * (n0 + n1 + n2 + n3 + n4);
	};

	ForgeNoise.prototype.generateTiling2D = function(x, y, options = {}) {
		let periodX = options.periodX || 256, periodY = options.periodY || 256;
		let xWrapped = x % periodX, yWrapped = y % periodY;
		return this.generate2D(xWrapped, yWrapped);
	};

	ForgeNoise.prototype.warp2D = function(x, y, options = {}) {
		let warpStrength = options.warpStrength || 1.0;
		let warpScale = options.warpScale || 1.0;
		let iterations = options.iterations || 1;
		let qx = x, qy = y;
		for (let i = 0; i < iterations; i++) {
			let offsetX = this.generate2D(qx * warpScale, qy * warpScale);
			let offsetY = this.generate2D(qx * warpScale + 5.2, qy * warpScale + 1.3);
			qx += offsetX * warpStrength;
			qy += offsetY * warpStrength;
		}
		return this.generate2D(qx, qy);
	};

	ForgeNoise.prototype.generate2D01 = function(x, y) { return (this.generate2D(x, y) + 1) * 0.5; };
	ForgeNoise.prototype.generate3D01 = function(x, y, z) { return (this.generate3D(x, y, z) + 1) * 0.5; };
	ForgeNoise.prototype.generateSimplex2D01 = function(x, y) { return (this.generateSimplex2D(x, y) + 1) * 0.5; };
	ForgeNoise.prototype.generateSimplex3D01 = function(x, y, z) { return (this.generateSimplex3D(x, y, z) + 1) * 0.5; };
	ForgeNoise.prototype.generateValue2D01 = function(x, y) { return this.generateValue2D(x, y); };
	ForgeNoise.prototype.generateValue3D01 = function(x, y, z) { return this.generateValue3D(x, y, z); };
	ForgeNoise.prototype.generateWorley2D01 = function(x, y) { return this.generateWorley2D(x, y) / 1.414; };
	ForgeNoise.prototype.generateWorley3D01 = function(x, y, z) { return this.generateWorley3D(x, y, z) / 1.732; };
	ForgeNoise.prototype.generate4D01 = function(x, y, z, w) { return (this.generate4D(x, y, z, w) + 1) * 0.5; };
	ForgeNoise.prototype.generateSimplex4D01 = function(x, y, z, w) { return (this.generateSimplex4D(x, y, z, w) + 1) * 0.5; };
	ForgeNoise.prototype.generateTiling2D01 = function(x, y, options) { return (this.generateTiling2D(x, y, options) + 1) * 0.5; };

	if (typeof module !== 'undefined' && module.exports) {
		module.exports = ForgeNoise;
	} else {
		root.ForgeNoise = ForgeNoise;
	}
})(this);
