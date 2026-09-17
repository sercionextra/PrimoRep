/**
 * =========================================================================================
 * CORE_MATH_ELIOSTATO_(v.4.3).js
 * 
 * Master Blueprint Ingegneristico - Libreria Standalone di Calcolo Fisico-Matematico
 * Robot Parallelo a Cavi (CDPR) Sottovincolato a 3 Cavi per Eliostato Solare
 * 
 * Release v.4.3:
 *  - Calcolo analitico dell'orario locale esatto (HH:MM) di Alba e Tramonto con Equazione del Tempo.
 *  - Gestione dinamica dei regimi di fuso orario: Auto (UE), Ora Solare (CET) e Ora Legale (CEST).
 *  - Navigatore polare completo con trasformazioni cartesiane per bersaglio a parete 360°.
 *  - Correzione trigonometrica esatta della taratura a due punti (R_z(-gamma)).
 * =========================================================================================
 */

(function (global, factory) {
    if (typeof exports === 'object' && typeof module !== 'undefined') {
        module.exports = factory();
    } else if (typeof define === 'function' && define.amd) {
        define(factory);
    } else {
        global.HeliostatCDPR = factory();
    }
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';

    // =========================================================================
    // 1. COSTANTI FISICHE E PARAMETRI GEOMETRICI NOMINALI
    // =========================================================================
    const DEFAULTS = {
        R_A: 0.30,                 // Raggio pulegge fisse superiori (m)
        r_b: 0.15,                 // Raggio vertici specchio mobile (m)
        PLATE_THICKNESS: 0.005,    // Spessore sandwich policarbonato + metallo (m)
        H_NOMINAL: 0.35,           // Quota baricentro nominale (m)
        H_MAX_DEPTH: 0.42,         // Quota massima Z-adattativa (m)
        M_TOTAL: 1.0,              // Massa complessiva specchio + piastra zavorra (kg)
        G_ACC: 9.80665,            // Accelerazione gravitazionale standard (m/s^2)
        MU_CAPSTAN: 0.13,          // Coefficiente attrito Dyneema PU su alluminio anodizzato
        WRAP_ANGLE: 6.0 * Math.PI, // 3 giri completi attorno al tamburo (rad)
        M_CW: 0.19,                // Massa contrappeso passivo per cavo (kg)
        R_CYLINDER: 0.0125,        // Raggio tamburo cabestano (m, diametro 25mm)
        STEPS_PER_REV: 3200,       // 200 passi/giro x 16 microstep
        TARGET_DISTANCE: 40.0,     // Distanza nominale bersaglio (m)
        CABLE_ANGLES: [            // Angoli geometrici pulegge (rad)
            Math.PI / 2.0,         // Cavo 1: Nord (90°)
            7.0 * Math.PI / 6.0,   // Cavo 2: Sud-Ovest (210°)
            11.0 * Math.PI / 6.0   // Cavo 3: Sud-Est (330°)
        ]
    };

    DEFAULTS.T_MIN_SLIP = (DEFAULTS.M_CW * DEFAULTS.G_ACC) / Math.exp(DEFAULTS.MU_CAPSTAN * DEFAULTS.WRAP_ANGLE);

    // =========================================================================
    // 2. LIVELLO 1: VMATH (ALGEBRA VETTORIALE/MATRICIALE ZERO-GC RIGOROSA)
    // =========================================================================
    const VMath = {
        create: (x = 0, y = 0, z = 0) => [x, y, z],
        clone: (v, out = [0, 0, 0]) => { out[0] = v[0]; out[1] = v[1]; out[2] = v[2]; return out; },
        set: (out, x, y, z) => { out[0] = x; out[1] = y; out[2] = z; return out; },

        add: (a, b, out = [0, 0, 0]) => {
            out[0] = a[0] + b[0]; out[1] = a[1] + b[1]; out[2] = a[2] + b[2];
            return out;
        },
        sub: (a, b, out = [0, 0, 0]) => {
            out[0] = a[0] - b[0]; out[1] = a[1] - b[1]; out[2] = a[2] - b[2];
            return out;
        },
        scale: (a, s, out = [0, 0, 0]) => {
            out[0] = a[0] * s; out[1] = a[1] * s; out[2] = a[2] * s;
            return out;
        },
        dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
        cross: (a, b, out = [0, 0, 0]) => {
            const x = a[1] * b[2] - a[2] * b[1];
            const y = a[2] * b[0] - a[0] * b[2];
            const z = a[0] * b[1] - a[1] * b[0];
            out[0] = x; out[1] = y; out[2] = z;
            return out;
        },
        norm: (a) => Math.hypot(a[0], a[1], a[2]),
        dist: (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]),
        normalize: (a, out = [0, 0, 0]) => {
            const n = Math.hypot(a[0], a[1], a[2]);
            if (n > 1e-12) {
                out[0] = a[0] / n; out[1] = a[1] / n; out[2] = a[2] / n;
            } else {
                out[0] = 0; out[1] = 0; out[2] = 1;
            }
            return out;
        },
        angleBetween: (a, b) => {
            const d = VMath.dot(a, b) / (Math.max(1e-12, VMath.norm(a) * VMath.norm(b)));
            return Math.acos(Math.max(-1.0, Math.min(1.0, d)));
        },

        rotationMatrixMollerHughes: (n_des, psi = 0) => {
            const n = VMath.normalize(n_des);
            const c = n[2];
            const v0 = -n[1], v1 = n[0];
            const factor = c > -0.99999 ? 1.0 / (1.0 + c) : 0.0;

            const R0 = [
                [1.0 - v1 * v1 * factor,  v0 * v1 * factor,        v1],
                [v0 * v1 * factor,        1.0 - v0 * v0 * factor, -v0],
                [-v1,                     v0,                      c ]
            ];

            const c_psi = Math.cos(psi);
            const s_psi = Math.sin(psi);

            const R = [ [0,0,0], [0,0,0], [0,0,0] ];
            for (let i = 0; i < 3; i++) {
                R[i][0] =  R0[i][0] * c_psi + R0[i][1] * s_psi;
                R[i][1] = -R0[i][0] * s_psi + R0[i][1] * c_psi;
                R[i][2] =  R0[i][2];
            }
            return R;
        },

        matVecMul3: (R, v, out = [0, 0, 0]) => {
            const x = R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2];
            const y = R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2];
            const z = R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2];
            out[0] = x; out[1] = y; out[2] = z;
            return out;
        },

        solveLinear3x3: (A, b, out = [0, 0, 0]) => {
            const d = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
                    - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
                    + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
            if (Math.abs(d) < 1e-12) return null;
            const invD = 1.0 / d;

            out[0] = (b[0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
                    - A[0][1]*(b[1]*A[2][2] - A[1][2]*b[2])
                    + A[0][2]*(b[1]*A[2][1] - A[1][1]*b[2])) * invD;

            out[1] = (A[0][0]*(b[1]*A[2][2] - A[1][2]*b[2])
                    - b[0]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
                    + A[0][2]*(A[1][0]*b[2] - b[1]*A[2][0])) * invD;

            out[2] = (A[0][0]*(A[1][1]*b[2] - b[1]*A[2][1])
                    - A[0][1]*(A[1][0]*b[2] - b[1]*A[2][0])
                    + b[0]*(A[1][0]*A[2][1] - A[1][1]*A[2][0])) * invD;

            return out;
        }
    };

    // =========================================================================
    // 3. LIVELLO 2: ASTRONOMIA SOLARE PSA CON CALCOLO ORARI ALBA/TRAMONTO (v.4.3)
    // =========================================================================
    const AstronomySolver = {
        dateToJulianDay: (date) => {
            let year = date.getUTCFullYear();
            let month = date.getUTCMonth() + 1;
            const day = date.getUTCDate() + (date.getUTCHours() + (date.getUTCMinutes() + (date.getUTCSeconds() + date.getUTCMilliseconds()/1000.0)/60.0)/60.0)/24.0;
            if (month <= 2) { year -= 1; month += 12; }
            const A = Math.floor(year / 100);
            const B = 2 - A + Math.floor(A / 4);
            return Math.floor(365.25 * (year + 4716)) + Math.floor(30.6001 * (month + 1)) + day + B - 1524.5;
        },

        calculateSunPosition: (dateOrEpoch, latDeg, lonDeg) => {
            const date = (dateOrEpoch instanceof Date) ? dateOrEpoch : new Date(dateOrEpoch);
            const JD = AstronomySolver.dateToJulianDay(date);
            const n = JD - 2451545.0;

            const deg2rad = Math.PI / 180.0;
            const rad2deg = 180.0 / Math.PI;

            let L = (280.460 + 0.9856474 * n) % 360.0;
            if (L < 0) L += 360.0;
            let g = (357.528 + 0.9856003 * n) % 360.0;
            if (g < 0) g += 360.0;

            const lambda = L + 1.915 * Math.sin(g * deg2rad) + 0.020 * Math.sin(2.0 * g * deg2rad);
            const epsilon = 23.439 - 0.0000004 * n;

            const sinDelta = Math.sin(epsilon * deg2rad) * Math.sin(lambda * deg2rad);
            const delta = Math.asin(Math.max(-1.0, Math.min(1.0, sinDelta)));
            const alpha = Math.atan2(Math.cos(epsilon * deg2rad) * Math.sin(lambda * deg2rad), Math.cos(lambda * deg2rad));

            let gmst = (280.46061837 + 360.98564736629 * n) % 360.0;
            if (gmst < 0) gmst += 360.0;
            
            let lmstDeg = (gmst + lonDeg) % 360.0;
            if (lmstDeg < 0) lmstDeg += 360.0;
            const lmst = lmstDeg * deg2rad;

            let omega = lmst - alpha;
            omega = Math.atan2(Math.sin(omega), Math.cos(omega));

            const phi = latDeg * deg2rad;
            const cosPhi = Math.cos(phi);
            const sinAlpha = Math.sin(phi) * Math.sin(delta) + cosPhi * Math.cos(delta) * Math.cos(omega);
            const elevation = Math.asin(Math.max(-1.0, Math.min(1.0, sinAlpha)));

            const cosElev = Math.cos(elevation);
            let azimuth = Math.PI;
            if (cosElev > 1e-6 && Math.abs(cosPhi) > 1e-6) {
                const cosGamma = (Math.sin(delta) - Math.sin(phi) * Math.sin(elevation)) / (cosPhi * cosElev);
                azimuth = Math.acos(Math.max(-1.0, Math.min(1.0, cosGamma)));
                if (Math.sin(omega) > 0) {
                    azimuth = 2.0 * Math.PI - azimuth;
                }
            }

            // Vettore Sole nel sistema geografico ENU [Est, Nord, Zenit]
            const sunVector = [
                Math.cos(elevation) * Math.sin(azimuth),
                Math.cos(elevation) * Math.cos(azimuth),
                Math.sin(elevation)
            ];

            return {
                elevationRad: elevation,
                azimuthRad: azimuth,
                elevationDeg: elevation * rad2deg,
                azimuthDeg: azimuth * rad2deg,
                sunVector: VMath.normalize(sunVector)
            };
        },

        // Calcolo esatto di Azimut e Ora Locale di Alba e Tramonto con Equazione del Tempo
        computeRiseSetData: (year, dayOfYear, latDeg, lonDeg, tzOffsetHours = 2.0) => {
            const midDayDate = new Date(Date.UTC(year, 0, dayOfYear, 12, 0, 0));
            const JD = AstronomySolver.dateToJulianDay(midDayDate);
            const n = JD - 2451545.0;
            const deg2rad = Math.PI / 180.0;

            let L = (280.460 + 0.9856474 * n) % 360.0;
            if (L < 0) L += 360.0;
            let g = (357.528 + 0.9856003 * n) % 360.0;
            if (g < 0) g += 360.0;
            const lambda = L + 1.915 * Math.sin(g * deg2rad) + 0.020 * Math.sin(2.0 * g * deg2rad);
            const epsilon = 23.439 - 0.0000004 * n;
            const delta = Math.asin(Math.sin(epsilon * deg2rad) * Math.sin(lambda * deg2rad));
            const alpha = Math.atan2(Math.cos(epsilon * deg2rad) * Math.sin(lambda * deg2rad), Math.cos(lambda * deg2rad));

            // Equazione del Tempo in minuti
            let alphaDeg = alpha * (180.0 / Math.PI);
            if (alphaDeg < 0) alphaDeg += 360.0;
            const eotMinutes = 4.0 * (L - alphaDeg);

            const phi = latDeg * deg2rad;
            const cosPhi = Math.cos(phi);
            const h0 = -0.8333 * deg2rad; // Rifrazione atmosferica + diametro semidisco solare
            
            let azRise = 90.0, azSet = 270.0;
            let timeRiseStr = "--:--", timeSetStr = "--:--";

            if (Math.abs(cosPhi) > 1e-6) {
                const cosH0 = (Math.sin(h0) - Math.sin(phi) * Math.sin(delta)) / (cosPhi * Math.cos(delta));
                if (cosH0 >= -1.0 && cosH0 <= 1.0) {
                    const H0_rad = Math.acos(cosH0);
                    const H0_hours = H0_rad * (12.0 / Math.PI);

                    const cosAz = (Math.sin(delta) - Math.sin(phi) * Math.sin(h0)) / (cosPhi * Math.cos(h0));
                    const azRiseRad = Math.acos(Math.max(-1.0, Math.min(1.0, cosAz)));
                    azRise = azRiseRad * (180.0 / Math.PI);
                    azSet = 360.0 - azRise;

                    // Transito solare (Mezzogiorno solare) in UTC
                    const transitUtcHours = 12.0 - (lonDeg / 15.0) - (eotMinutes / 60.0);
                    const riseUtc = transitUtcHours - H0_hours;
                    const setUtc = transitUtcHours + H0_hours;

                    const toHM = (decHours) => {
                        let localH = (decHours + tzOffsetHours + 24.0) % 24.0;
                        let h = Math.floor(localH);
                        let m = Math.floor((localH % 1.0) * 60.0);
                        return `${String(h).padStart(2,'0')}:${String(m).padStart(2,'0')}`;
                    };

                    timeRiseStr = toHM(riseUtc);
                    timeSetStr = toHM(setUtc);
                }
            }
            return { azRise, azSet, timeRiseStr, timeSetStr, deltaRad: delta };
        }
    };

    // =========================================================================
    // 4. LIVELLO 3: CINEMATICA CDPR (MODELLO A & IK Z-ADATTATIVA)
    // =========================================================================
    class CDPRKinematics {
        constructor(config = {}) {
            this.cfg = Object.assign({}, DEFAULTS, config);
            this.L_SIDE = this.cfg.r_b * Math.sqrt(3.0);
            this.L_0_NOMINAL = Math.sqrt(Math.pow(this.cfg.R_A - this.cfg.r_b, 2) + Math.pow(this.cfg.H_NOMINAL, 2));

            this.A_pts = this.cfg.CABLE_ANGLES.map(ang => [
                this.cfg.R_A * Math.cos(ang),
                this.cfg.R_A * Math.sin(ang),
                0.0
            ]);

            this.b_pts = this.cfg.CABLE_ANGLES.map(ang => [
                this.cfg.r_b * Math.cos(ang),
                this.cfg.r_b * Math.sin(ang),
                0.0
            ]);

            this._ikScratch = {
                r: [[0,0,0], [0,0,0], [0,0,0]],
                p: [[0,0,0], [0,0,0], [0,0,0]],
                u: [[0,0,0], [0,0,0], [0,0,0]],
                L_vec: [0,0,0],
                tmp: [0,0,0],
                J: [ [0,0,0], [0,0,0], [0,0,0] ],
                A: [ [0,0,0], [0,0,0], [0,0,0] ],
                g: [0, 0, 0],
                dx: [0, 0, 0],
                statePert: [0, 0, 0],
                trialState: [0, 0, 0]
            };
        }

        computeStaticTensions(u) {
            const u2_x_u3 = VMath.cross(u[1], u[2]);
            const detU = VMath.dot(u[0], u2_x_u3);
            if (Math.abs(detU) < 1e-6) return { T: [0, 0, 0], detU: 0, isValid: false, isGripOk: false };

            const W = this.cfg.M_TOTAL * this.cfg.G_ACC;
            const T = [
                W * u2_x_u3[2] / detU,
                W * VMath.cross(u[2], u[0])[2] / detU,
                W * VMath.cross(u[0], u[1])[2] / detU
            ];

            const isValid = T[0] > 0.05 && T[1] > 0.05 && T[2] > 0.05;
            const isGripOk = T[0] >= this.cfg.T_MIN_SLIP && T[1] >= this.cfg.T_MIN_SLIP && T[2] >= this.cfg.T_MIN_SLIP;
            return { T, detU, isValid, isGripOk };
        }

        solveInverseKinematics(n_des, zTarget = -this.cfg.H_NOMINAL) {
            const sc = this._ikScratch;
            const nd = VMath.normalize(n_des);

            const tiltAngleRad = Math.acos(Math.max(-1.0, Math.min(1.0, nd[2])));
            const tiltDeg = tiltAngleRad * (180.0 / Math.PI);

            let effectiveZ = zTarget;
            if (tiltDeg > 28.0) {
                const normTilt = Math.min(1.0, (tiltDeg - 28.0) / 22.0);
                effectiveZ = zTarget - (normTilt * normTilt) * (this.cfg.H_MAX_DEPTH - Math.abs(zTarget));
            }

            const evalResiduals = (state, curZ) => {
                const xG = state[0], yG = state[1], psi = state[2];
                const R = VMath.rotationMatrixMollerHughes(nd, psi);

                for (let i = 0; i < 3; i++) {
                    VMath.matVecMul3(R, this.b_pts[i], sc.r[i]);
                    sc.p[i][0] = xG + sc.r[i][0];
                    sc.p[i][1] = yG + sc.r[i][1];
                    sc.p[i][2] = curZ + sc.r[i][2];

                    VMath.sub(this.A_pts[i], sc.p[i], sc.L_vec);
                    VMath.normalize(sc.L_vec, sc.u[i]);
                }

                const resT = this.computeStaticTensions(sc.u);
                
                let tauX = 0.0, tauY = 0.0, tauZ = 0.0;
                for (let i = 0; i < 3; i++) {
                    VMath.cross(sc.r[i], sc.u[i], sc.tmp);
                    tauX += resT.T[i] * sc.tmp[0];
                    tauY += resT.T[i] * sc.tmp[1];
                    tauZ += resT.T[i] * sc.tmp[2];
                }

                return {
                    F: [tauX, tauY, tauZ],
                    T: resT.T,
                    isValid: resT.isValid,
                    isGripOk: resT.isGripOk,
                    detU: resT.detU,
                    R
                };
            };

            let state = [0.065 * nd[0], 0.065 * nd[1], 0.0];
            let lambda = 1e-3;
            const EPS = 1e-4;
            let converged = false;
            let finalIter = 0;

            for (let iter = 0; iter < 16; iter++) {
                finalIter = iter;
                const res0 = evalResiduals(state, effectiveZ);
                const normF = Math.hypot(res0.F[0], res0.F[1], res0.F[2]);

                if (normF < 1e-5 && res0.isValid) {
                    converged = true;
                    break;
                }

                for (let j = 0; j < 3; j++) {
                    sc.statePert[0] = state[0];
                    sc.statePert[1] = state[1];
                    sc.statePert[2] = state[2];
                    sc.statePert[j] += EPS;
                    const resPert = evalResiduals(sc.statePert, effectiveZ);
                    for (let i = 0; i < 3; i++) {
                        sc.J[i][j] = (resPert.F[i] - res0.F[i]) / EPS;
                    }
                }

                for (let i = 0; i < 3; i++) {
                    for (let j = 0; j < 3; j++) {
                        let dot = 0.0;
                        for (let k = 0; k < 3; k++) dot += sc.J[k][i] * sc.J[k][j];
                        sc.A[i][j] = dot;
                    }
                    sc.A[i][i] += lambda * (sc.A[i][i] > 1e-4 ? sc.A[i][i] : 1.0);

                    let dotG = 0.0;
                    for (let k = 0; k < 3; k++) dotG += sc.J[k][i] * res0.F[k];
                    sc.g[i] = dotG;
                }

                const dx = VMath.solveLinear3x3(sc.A, sc.g, sc.dx);
                if (!dx) {
                    lambda *= 5.0;
                    continue;
                }

                let stepSize = 1.0;
                let accepted = false;

                for (let b = 0; b < 4; b++) {
                    sc.trialState[0] = state[0] - dx[0] * stepSize;
                    sc.trialState[1] = state[1] - dx[1] * stepSize;
                    sc.trialState[2] = state[2] - dx[2] * stepSize;

                    const trialRes = evalResiduals(sc.trialState, effectiveZ);
                    const trialNormF = Math.hypot(trialRes.F[0], trialRes.F[1], trialRes.F[2]);

                    if (trialNormF < normF) {
                        state[0] = sc.trialState[0];
                        state[1] = sc.trialState[1];
                        state[2] = sc.trialState[2];
                        lambda = Math.max(1e-5, lambda * 0.5);
                        accepted = true;
                        break;
                    }
                    stepSize *= 0.5;
                }

                if (!accepted) lambda *= 4.0;
            }

            const finalEval = evalResiduals(state, effectiveZ);
            const lengths = [
                VMath.dist(this.A_pts[0], sc.p[0]),
                VMath.dist(this.A_pts[1], sc.p[1]),
                VMath.dist(this.A_pts[2], sc.p[2])
            ];

            return {
                success: converged || finalEval.isValid,
                pose: { x_G: state[0], y_G: state[1], z_G: effectiveZ, psi_rad: state[2] },
                lengths,
                vertices: [ [...sc.p[0]], [...sc.p[1]], [...sc.p[2]] ],
                tensions: finalEval.T,
                isGripOk: finalEval.isGripOk,
                iterations: finalIter,
                isZAdapted: Math.abs(effectiveZ - zTarget) > 1e-4
            };
        }

        relaxForwardKinematics(currentVertices, targetLengths, subSteps = 16) {
            const v = [ [...currentVertices[0]], [...currentVertices[1]], [...currentVertices[2]] ];
            const toP = [0,0,0], pG = [0,0,0], e1 = [0,0,0], e2 = [0,0,0], n_plane = [0,0,0], vY = [0,0,0], vX = [0,0,0];

            for (let i = 0; i < 3; i++) v[i][2] -= 0.0015;

            for (let s = 0; s < subSteps; s++) {
                for (let i = 0; i < 3; i++) {
                    VMath.sub(v[i], this.A_pts[i], toP);
                    const d = VMath.norm(toP);
                    const L = targetLengths[i];
                    if (d > L) {
                        const scale = L / d;
                        v[i][0] = this.A_pts[i][0] + toP[0] * scale;
                        v[i][1] = this.A_pts[i][1] + toP[1] * scale;
                        v[i][2] = this.A_pts[i][2] + toP[2] * scale;
                    }
                }

                pG[0] = (v[0][0] + v[1][0] + v[2][0]) / 3.0;
                pG[1] = (v[0][1] + v[1][1] + v[2][1]) / 3.0;
                pG[2] = (v[0][2] + v[1][2] + v[2][2]) / 3.0;

                VMath.sub(v[1], v[0], e1);
                VMath.sub(v[2], v[0], e2);
                VMath.cross(e1, e2, n_plane);
                if (n_plane[2] < 0) VMath.scale(n_plane, -1.0, n_plane);
                VMath.normalize(n_plane, n_plane);

                VMath.sub(v[0], pG, vY);
                const proj = VMath.dot(vY, n_plane);
                vY[0] -= n_plane[0] * proj;
                vY[1] -= n_plane[1] * proj;
                vY[2] -= n_plane[2] * proj;
                VMath.normalize(vY, vY);
                VMath.cross(vY, n_plane, vX);

                for (let i = 0; i < 3; i++) {
                    const bx = this.b_pts[i][0], by = this.b_pts[i][1];
                    v[i][0] = pG[0] + vX[0] * bx + vY[0] * by;
                    v[i][1] = pG[1] + vX[1] * bx + vY[1] * by;
                    v[i][2] = pG[2] + vX[2] * bx + vY[2] * by;
                }
            }

            VMath.sub(v[1], v[0], e1);
            VMath.sub(v[2], v[0], e2);
            const n_eff = VMath.cross(e1, e2);
            if (n_eff[2] < 0) VMath.scale(n_eff, -1.0, n_eff);
            VMath.normalize(n_eff, n_eff);

            return { vertices: v, pG, n_eff };
        }
    }

    // =========================================================================
    // 5. LIVELLO 4: TARATURA OTTICA A 2 PUNTI RADDRIZZATA
    // =========================================================================
    const CalibrationEngine = {
        invertTargetVector: (s, n) => {
            const sNorm = VMath.normalize(s);
            const nNorm = VMath.normalize(n);
            const dot = VMath.dot(sNorm, nNorm);
            return VMath.normalize(VMath.sub(VMath.scale(nNorm, 2.0 * dot), sNorm));
        },

        solveTwoPointCalibration: (p1, p2, latDeg, lonDeg) => {
            const s1_geo = AstronomySolver.calculateSunPosition(p1.epoch, latDeg, lonDeg).sunVector;
            const s2_geo = AstronomySolver.calculateSunPosition(p2.epoch, latDeg, lonDeg).sunVector;

            let bestGamma = 0.0;
            let minResidual = 1e9;

            const COARSE_STEP = (1.0 * Math.PI) / 180.0;
            for (let gamma = 0; gamma < 2.0 * Math.PI; gamma += COARSE_STEP) {
                const c = Math.cos(gamma);
                const s = Math.sin(gamma);
                const s1_b = [ s1_geo[0]*c - s1_geo[1]*s,  s1_geo[0]*s + s1_geo[1]*c, s1_geo[2] ];
                const s2_b = [ s2_geo[0]*c - s2_geo[1]*s,  s2_geo[0]*s + s2_geo[1]*c, s2_geo[2] ];

                const t1 = CalibrationEngine.invertTargetVector(s1_b, p1.n_eff);
                const t2 = CalibrationEngine.invertTargetVector(s2_b, p2.n_eff);

                let err = 1.0 - VMath.dot(t1, t2);
                if (t1[2] < -0.05 || t2[2] < -0.05) err += 10.0;

                if (err < minResidual) {
                    minResidual = err;
                    bestGamma = gamma;
                }
            }

            const FINE_STEP = (0.01 * Math.PI) / 180.0;
            const gStart = bestGamma - (1.5 * Math.PI) / 180.0;
            const gEnd   = bestGamma + (1.5 * Math.PI) / 180.0;

            for (let gamma = gStart; gamma <= gEnd; gamma += FINE_STEP) {
                const c = Math.cos(gamma);
                const s = Math.sin(gamma);
                const s1_b = [ s1_geo[0]*c - s1_geo[1]*s,  s1_geo[0]*s + s1_geo[1]*c, s1_geo[2] ];
                const s2_b = [ s2_geo[0]*c - s2_geo[1]*s,  s2_geo[0]*s + s2_geo[1]*c, s2_geo[2] ];

                const t1 = CalibrationEngine.invertTargetVector(s1_b, p1.n_eff);
                const t2 = CalibrationEngine.invertTargetVector(s2_b, p2.n_eff);

                let err = 1.0 - VMath.dot(t1, t2);
                if (t1[2] < -0.05 || t2[2] < -0.05) err += 10.0;

                if (err < minResidual) {
                    minResidual = err;
                    bestGamma = gamma;
                }
            }

            const c_opt = Math.cos(bestGamma);
            const s_opt = Math.sin(bestGamma);
            const s1_opt = [ s1_geo[0]*c_opt - s1_geo[1]*s_opt, s1_geo[0]*s_opt + s1_geo[1]*c_opt, s1_geo[2] ];
            const s2_opt = [ s2_geo[0]*c_opt - s2_geo[1]*s_opt, s2_geo[0]*s_opt + s2_geo[1]*c_opt, s2_geo[2] ];

            const t1_opt = CalibrationEngine.invertTargetVector(s1_opt, p1.n_eff);
            const t2_opt = CalibrationEngine.invertTargetVector(s2_opt, p2.n_eff);
            const t_target = VMath.normalize(VMath.add(t1_opt, t2_opt));

            let finalDeg = (bestGamma * 180.0 / Math.PI) % 360.0;
            if (finalDeg < 0) finalDeg += 360.0;

            return {
                gammaBaseRad: bestGamma,
                gammaBaseDeg: finalDeg,
                targetVectorBase: t_target,
                residualError: minResidual
            };
        }
    };

    // =========================================================================
    // 6. LIVELLO 5: GUI HELPERS & ADATTATORE THREE.JS
    // =========================================================================
    const GUIHelpers = {
        computeThreeOrthonormalBasis: (vertices, pG, n_eff) => {
            const vY = VMath.normalize(VMath.sub(vertices[0], pG));
            const vZ = VMath.normalize(n_eff);
            const vX = VMath.normalize(VMath.cross(vY, vZ));
            const vY_reortho = VMath.normalize(VMath.cross(vZ, vX));
            return { vX, vY: vY_reortho, vZ };
        },

        computeMirrorProjectionTrace: (pG, n_des, n_eff, r_b = DEFAULTS.r_b) => {
            const nd = VMath.normalize(n_des);
            const ne = VMath.normalize(n_eff);
            const dot = VMath.dot(nd, ne);

            const vProj = VMath.sub(nd, VMath.scale(ne, dot));
            const sinDeltaTheta = VMath.norm(vProj);

            const pStart = VMath.add(pG, VMath.scale(ne, 0.002));
            const displayLength = sinDeltaTheta * r_b;
            const pEnd = VMath.clone(pStart);

            if (sinDeltaTheta > 1e-5) {
                const dir = VMath.scale(vProj, 1.0 / sinDeltaTheta);
                pEnd[0] += dir[0] * displayLength;
                pEnd[1] += dir[1] * displayLength;
                pEnd[2] += dir[2] * displayLength;
            }

            return {
                pStart,
                pEnd,
                lengthMeters: displayLength,
                lengthMm: displayLength * 1000.0,
                vProj
            };
        },

        formatTelemetryData: (n_eff, n_des, cableLengths) => {
            const angleRad = VMath.angleBetween(n_eff, n_des);
            const errDeg = angleRad * (180.0 / Math.PI);

            let status = 'OFF TARGET';
            let badgeClass = 'error-bad';

            if (errDeg < 0.25) {
                status = 'LOCK (ALLINEATO)';
                badgeClass = 'error-good';
            } else if (errDeg < 2.5) {
                status = 'APPROACH';
                badgeClass = 'error-mid';
            }

            return {
                errorDeg: errDeg,
                errorDegFormatted: errDeg.toFixed(2) + '°',
                status,
                badgeClass,
                nEffFormatted: `[${n_eff[0].toFixed(2)}, ${n_eff[1].toFixed(2)}, ${n_eff[2].toFixed(2)}]`,
                nDesFormatted: `[${n_des[0].toFixed(2)}, ${n_des[1].toFixed(2)}, ${n_des[2].toFixed(2)}]`,
                cableLengthsMm: cableLengths.map(l => (l * 1000.0).toFixed(1))
            };
        }
    };

    // =========================================================================
    // 7. MASTER FACADE ENGINE
    // =========================================================================
    class HeliostatMasterEngine {
        constructor(config = {}) {
            this.config = Object.assign({}, DEFAULTS, config);
            this.kinematics = new CDPRKinematics(this.config);
            this.astronomy = AstronomySolver;
            this.calibration = CalibrationEngine;
            this.helpers = GUIHelpers;
            this.VMath = VMath;
            this.version = "4.3.0";
        }
    }

    return {
        VMath,
        AstronomySolver,
        CDPRKinematics,
        CalibrationEngine,
        GUIHelpers,
        HeliostatMasterEngine,
        VERSION: "4.3.0"
    };
}));