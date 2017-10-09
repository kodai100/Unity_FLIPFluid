using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Kodai.FLIP.CPU {

    public class FLIP : MonoBehaviour {

        public enum Cell {
            EMPTY, FLUID, SOLID
        }

        public int n = 32;
        public float dt = 1 / 120f;
        public float mass = 1f;
        public Vector2 iv = new Vector2(0.1f, 0);
        public Vector2 force = new Vector2(0, -9.8f);

        int y;

        private float[] u;
        private float[] ux;
        private float[] m;
        private float[] p;
        private float[] divu;

        private float[] px, ipx;
        private float[] py, ipy;
        private float[] pvx;
        private float[] pvy;

        private Cell[] flags;

        private int numParticles;


        void Start() {
            Initialize();
        }


        void Update() {
            particles2grid();

            for (int k = 0; k < n * (n + 1); ++k) {
                if (m[y + k] > 1e-8) u[y + k] += force.y * dt * (n - 2);
                if (m[k] > 1e-8) u[k] += force.x * dt * (n - 2);
            }

            pressureSolve();
            updateAndAdvect(0.96f);
        }

        private void OnDrawGizmos() {
            if (!Application.isPlaying) return;
            for (int i = 0; i < px.Length; i++) {
                Gizmos.DrawWireSphere(new Vector3(px[i], py[i], 0), 0.1f);
            }
        }

        void Initialize() {
            y = n * (n + 1);
            u = new float[2 * n * (n + 1)];
            ux = new float[2 * n * (n + 1)];
            m = new float[2 * n * (n + 1)];
            p = new float[n * n];
            divu = new float[n * n];
            flags = new Cell[n * n];


            numParticles = (n - 2) / 4 * (n - 4) * 4;
            Debug.Log(numParticles);

            px = new float[numParticles];
            py = new float[numParticles];
            ipx = new float[numParticles];
            ipy = new float[numParticles];
            pvx = new float[numParticles];
            pvy = new float[numParticles];

            for (int i = 0; i < numParticles; i++)
            {
                pvx[i] = iv.x;
                pvy[i] = iv.y;
            }

            int idx = 0;
            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    if (idx < numParticles && i - 1 < (n - 2) / 4) {
                        for (int jj = 0; jj < 2; ++jj) {
                            for (int ii = 0; ii < 2; ++ii) {
                                px[idx] = i + 0.25f + ii * 0.5f;
                                ipx[idx] = px[idx];
                                py[idx++] = j + 0.25f + jj * 0.5f;
                                ipy[idx-1] = py[idx-1];
                            }
                        }
                    }
                }
            }
        }

        float InterpolateX(int i, int j, float px, float py)
        {
            return (1f - Mathf.Abs(i - px)) * (1f - Mathf.Abs(j + 0.5f - py));
        }

        float InterpolateY(int i, int j, float px, float py)
        {
            return (1f - Mathf.Abs(i + 0.5f - px)) * (1f - Mathf.Abs(j - py));
        }

    void particles2grid() {

            // memset u 0
            // memset m 0
            for(int i = 0; i < u.Length; i++) {
                u[i] = 0;
                m[i] = 0;
            }

            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    flags[i + j * n] = (i == 0 || j == 0 || i == n - 1 || j == n - 1) ? Cell.SOLID : Cell.EMPTY;
                }
            }
            

            for (int k = 0; k < numParticles; ++k) {
                int i = (int)px[k], j = (int)py[k], fi = (int)(px[k] - 0.5f), fj = (int)(py[k] - 0.5f);
                flags[i + j * n] = Cell.FLUID;
                for (int jj = fj; jj < fj + 2; ++jj) {
                    for (int ii = i; ii < i + 2; ++ii) {
                        u[IDXX(ii, jj)] += pvx[k] * InterpolateX(ii, jj, px[k], py[k]);
                        m[IDXX(ii, jj)] += mass * InterpolateX(ii, jj, px[k], py[k]);
                    }
                }
                for (int jj = j; jj < j + 2; ++jj) {
                    for (int ii = fi; ii < fi + 2; ++ii) {
                        u[IDXY(ii, jj)] += pvy[k] * InterpolateY(ii, jj, px[k], py[k]);
                        m[IDXY(ii, jj)] += mass * InterpolateY(ii, jj, px[k], py[k]);
                    }
                }
            }

            for (int k = 0; k < 2 * n * (n + 1); ++k) {
                if (m[k] > 1e-8) {
                    m[k] = 1 / m[k]; u[k] *= m[k];
                }
            }

            // memcpy ux, u
            for (int i = 0; i < u.Length; i++) {
                ux[i] = u[i];
            }
        }

        void pressureSolve() {
            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    divu[i + j * n] = -(u[i + 1 + j * (n+1)] - u[i + j * (n+1)] + u[y + i + (j + 1) * n] - u[y + i + j * n]);
                }
            }

            for (int k = 1; k < n - 1; ++k) {
                divu[1 + k * n] -= u[1 + k * (n+1)];
                divu[n - 2 + k * n] += u[n - 1 + k * (n+1)];
                divu[k + n] -= u[y + k + n];
                divu[k + (n - 2) * n] += u[y + k + (n - 1) * n];
            }

            // memset(p, 0, sizeof(p));
            for(int i = 0; i < p.Length; i++) {
                p[i] = 0;
            }

            for (int k = 0; k < 10 * n * n; ++k) {    // Gauss-Seidel solve for pressure
                for (int j = 1; j < n - 1; ++j) {
                    for (int i = 1; i < n - 1; ++i) {
                        if ((int)flags[i + j * n] != 1) continue;
                        float l = i > 1 ? -1 : 0, r = i < n - 2 ? -1 : 0, b = j > 1 ? -1 : 0, t = j < n - 2 ? -1 : 0, c = -(l + r + b + t);
                        p[i + j * n] = 1 / c * (divu[i + j * n] - l * p[i - 1 + j * n] - r * p[i + 1 + j * n] - b * p[i + j * n - n] - t * p[i + j * n + n]);
                    }
                }
            }

            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n; ++i) {
                    u[i + j * (n + 1)] -= p[i + j * n] - p[i - 1 + j * n];
                }
            }
            for (int j = 1; j < n; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    u[y + i + j * n] -= p[i + j * n] - p[i + j * n - n];
                }
            }

            for (int k = 1; k < n - 1; ++k) {
                u[1 + k * (n + 1)] = u[n - 1 + k * (n+1)] = u[y + k + n] = u[y + k + (n - 1) * n] = 0;
            }
        }

        void updateAndAdvect(float flip) {
            for (int k = 0; k < numParticles; ++k) {
                int i = (int)(px[k]), j = (int)(py[k]), fi = (int)(px[k] - 0.5f), fj = (int)(py[k] - 0.5f);
                float vx = 0, vy = 0, vxOld = 0, vyOld = 0;

                for (int jj = fj; jj < fj + 2; ++jj) {
                    for (int ii = i; ii < i + 2; ++ii) {
                        vx += u[IDXX(ii, jj)] * InterpolateX(ii, jj, px[k], py[k]);
                        vxOld += ux[IDXX(ii, jj)] * InterpolateX(ii, jj, px[k], py[k]);
                    }
                }

                for (int jj = j; jj < j + 2; ++jj) {
                    for (int ii = fi; ii < fi + 2; ++ii) {
                        vy += u[IDXY(ii, jj)] * InterpolateY(ii, jj, px[k], py[k]);
                        vyOld += ux[IDXY(ii, jj)] * InterpolateY(ii, jj, px[k], py[k]);
                    }
                }

                pvx[k] = (1 - flip) * vx + flip * (pvx[k] + vx - vxOld);
                pvy[k] = (1 - flip) * vy + flip * (pvy[k] + vy - vyOld);

                // px[k] = Mathf.Min(Mathf.Max(px[k] + vx * dt, 1.001f), n - 1.001f);  // 境界
                // py[k] = Mathf.Min(Mathf.Max(py[k] + vy * dt, 1.001f), n - 1.001f);  // 境界

                if(px[k] + vx * dt < 1.001f)
                {
                    ReInitializeParticle(k);
                } else if(px[k] + vx * dt > n - 1.001f)
                {
                    ReInitializeParticle(k);
                }
                else
                {
                    px[k] = px[k] + vx * dt;
                }

                if (py[k] + vy * dt < 1.001f)
                {
                    ReInitializeParticle(k);
                }
                else if (py[k] + vy * dt > n - 1.001f)
                {
                    ReInitializeParticle(k);
                }
                else
                {
                    py[k] = py[k] + vy * dt;
                }
            }
        }

        void ReInitializeParticle(int k)
        {
            px[k] = 1.001f;
            py[k] = ipy[k];
            pvx[k] = iv.x;
            pvy[k] = iv.y;
        }

        int IDXX(int i, int j) {
            return (i + j * (n + 1));
        }

        int IDXY(int i, int j) {
            return y + i + j * n;
        }

    }

}