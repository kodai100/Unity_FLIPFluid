using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Kodai.FLIP.CPU {

    public enum Cell {
        EMPTY = 0, FLUID = 1, SOLID = 2
    }

    struct FluidParticle {
        public Vector2 pos;
        public Vector2 vel;
    }

    public class FLIP : MonoBehaviour {
        
        public int n = 32;
        public float dt = 1 / 120f;
        public float flip = 0.96f;

        public Color gizmoColor;
        public Material renderMat;

        private int y;

        FluidParticle[] particles;
        private float[] gridVel;
        private float[] gridVelSaved;
        private float[] gridMass;
        private float[] gridPressure;
        private float[] gridVelDivergence;
        private int[]  gridFlag;

        [SerializeField] private int numParticles;

        #region Renderer
        private Mesh mesh;
        private Vector3[] v;
        private int[] vIdx;
        #endregion
        
        #region MonoFuncs
        void Start() {
            InitializeGrid();
            InitializeParticle();
            InitializeRenderer();
        }

        
        void Update() {
            
            Transfer();
            AddForce();
            SolvePressure();
            Advect();

            Render();
            
        }

        private void OnDrawGizmos() {

            Gizmos.color = gizmoColor;

            Vector2 a = new Vector2(n, n);
            Gizmos.DrawWireCube(a/2, a);

            for (int i = 0; i < n; i++) {
                Gizmos.DrawLine(new Vector3(i, 0, 0), new Vector3(i, n, 0));
            }

            for (int i = 0; i < n; i++) {
                Gizmos.DrawLine(new Vector3(0, i, 0), new Vector3(n, i, 0));
            }
        }
        #endregion MonoFuncs

        #region PrivateFuncs

        #region Initializer
        void InitializeGrid() {
            // staggered grid
            y = n * (n + 1);
            gridVel = new float[2 * n * (n + 1)];
            gridVelSaved = new float[2 * n * (n + 1)];
            gridMass = new float[2 * n * (n + 1)];

            gridPressure = new float[n * n];
            gridVelDivergence = new float[n * n];
            gridFlag = new int[n * n];
        }

        void InitializeParticle() {
            numParticles = (n - 2) / 4 * (n - 4) * 4;
            particles = new FluidParticle[numParticles];

            int idx = 0;
            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    if (idx < numParticles && i - 1 < (n - 2) / 4) {
                        for (int jj = 0; jj < 2; ++jj) {
                            for (int ii = 0; ii < 2; ++ii) {
                                particles[idx].pos.x = i + 0.25f + ii * 0.5f;
                                particles[idx++].pos.y = j + 0.25f + jj * 0.5f;
                            }
                        }
                    }
                }
            }
        }

        void InitializeRenderer() {
            mesh = new Mesh();
            v = new Vector3[numParticles];
            vIdx = new int[numParticles];
        }
        #endregion Initializer

        #region GridProcess
        /// <summary>
        /// Reset grid quantities
        /// </summary>
        void ResetGrid() {
            // 2 * n * (n + 1)
            for (int i = 0; i < gridVel.Length; i++) {
                gridVel[i] = 0;
                gridMass[i] = 0;
            }

            // n * n
            for(int i = 0; i < gridPressure.Length; i++) {
                gridPressure[i] = 0;
                gridVelDivergence[i] = 0;
            }

            // 端の場合SOLID, それ以外の場合EMPTY
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    gridFlag[i + j * n] = (i == 0 || j == 0 || i == n - 1 || j == n - 1) ? (int)Cell.SOLID : (int)Cell.EMPTY;
                }
            }
        }

        /// <summary>
        /// Transfer particle velocities and masses to grid
        /// </summary>
        void Transfer() {

            ResetGrid();
            
            // Each particle
            for (int k = 0; k < numParticles; ++k) {
                // TODO 格子の大きさを可変に
                int i  = (int) particles[k].pos.x, 
                    j  = (int) particles[k].pos.y, 
                    fi = (int) (particles[k].pos.x - 0.5f), // 半分左にずらす
                    fj = (int) (particles[k].pos.y - 0.5f); // 半分下にずらす

                gridFlag[i + j * n] = (int)Cell.FLUID;  // 粒子の存在する格子は流体フラグ

                // X - velocity
                for (int jj = fj; jj < fj + 2; ++jj) {
                    for (int ii = i; ii < i + 2; ++ii) {
                        gridVel[CalcIndexX(ii, jj)] += particles[k].vel.x * InterpolateX(ii, jj, particles[k].pos);
                        gridMass[CalcIndexX(ii, jj)] += 1;// InterpolateX(ii, jj, particles[k].pos);
                    }
                }

                // Y - velocity
                for (int jj = j; jj < j + 2; ++jj) {
                    for (int ii = fi; ii < fi + 2; ++ii) {
                        gridVel[CalcIndexY(ii, jj)] += particles[k].vel.y * InterpolateY(ii, jj, particles[k].pos);
                        gridMass[CalcIndexY(ii, jj)] += InterpolateY(ii, jj, particles[k].pos);
                    }
                }
            }

            DebugConsole(gridMass);
            DebugConsole(gridFlag);
            

            for (int k = 0; k < 2 * n * (n + 1); ++k) {
                // 格子に質量があれば
                if (gridMass[k] > 1e-8) {
                    gridMass[k] = 1 / gridMass[k];   // 格子の質量を逆数にして
                    gridVel[k] *= gridMass[k];       // 速度に掛ける
                }
            }
            

            // Save current grid velocities
            for (int i = 0; i < gridVel.Length; i++) {
                gridVelSaved[i] = gridVel[i];
            }

        }

        void DebugConsole<T>(T[] array) {
            string str = "CPU: ";
            foreach (T tmp in array) {
                str += tmp + ", ";
            }
            Debug.Log(str);
        }

        /// <summary>
        /// Add acceleration(forces) to grid velocitiy
        /// </summary>
        void AddForce() {
            // y軸方向に重力を作用させる
            for (int k = 0; k < n * (n + 1); ++k) {
                if (gridMass[y + k] > 1e-8) gridVel[y + k] += -9.8f * dt * (n - 2);  // (n - 2) for simulation scale
            }
        }

        void SolvePressure() {
            
            // Calc divergence of grid velocities
            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    // 流入と流出の差の和のマイナス
                    gridVelDivergence[i + j * n] = -(gridVel[i + 1 + j * (n+1)] - gridVel[i + j * (n+1)] + gridVel[y + i + (j + 1) * n] - gridVel[y + i + j * n]);
                }
            }

            for (int k = 1; k < n - 1; ++k) {
                gridVelDivergence[1 + k * n] -= gridVel[1 + k * (n+1)];             // 左壁
                gridVelDivergence[n - 2 + k * n] += gridVel[n - 1 + k * (n+1)];     // 右壁
                gridVelDivergence[k + n] -= gridVel[y + k + n];                     // 下壁
                gridVelDivergence[k + (n - 2) * n] += gridVel[y + k + (n - 1) * n]; // 上壁
            }
            

            // Solve laplace ewuation using gauss-seidel
            for (int k = 0; k < 10 * n * n; ++k) {
                for (int j = 1; j < n - 1; ++j) {
                    for (int i = 1; i < n - 1; ++i) {
                        
                        if (gridFlag[i + j * n] != (int)Cell.FLUID) continue;

                        // 上下左右を判定し、枠だった場合は0
                        float l = i > 1 ? -1 : 0,       // left
                              r = i < n - 2 ? -1 : 0,   // right
                              b = j > 1 ? -1 : 0,       // bottom
                              t = j < n - 2 ? -1 : 0,   // top
                              c = -(l + r + b + t);     // center

                        // 圧力の計算
                        gridPressure[i + j * n] = 1 / c * (gridVelDivergence[i + j * n] - l * gridPressure[i - 1 + j * n] - r * gridPressure[i + 1 + j * n] - b * gridPressure[i + j * n - n] - t * gridPressure[i + j * n + n]);
                    }
                }
            }


            // 圧力差による速度の更新
            // x
            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n; ++i) {
                    gridVel[i + j * (n + 1)] -= gridPressure[i + j * n] - gridPressure[i - 1 + j * n];
                }
            }

            // y
            for (int j = 1; j < n; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    gridVel[y + i + j * n] -= gridPressure[i + j * n] - gridPressure[i + j * n - n];
                }
            }

            // 壁は0
            for (int k = 1; k < n - 1; ++k) {
                gridVel[1 + k * (n + 1)] = gridVel[n - 1 + k * (n+1)] = gridVel[y + k + n] = gridVel[y + k + (n - 1) * n] = 0;
            }
        }
        #endregion GridProcess

        void Advect() {

            // 各粒子について
            for (int k = 0; k < numParticles; ++k) {

                // 自身の格子インデックス
                int i  = (int) (particles[k].pos.x),
                    j  = (int) (particles[k].pos.y), 
                    fi = (int) (particles[k].pos.x - 0.5f), 
                    fj = (int) (particles[k].pos.y - 0.5f);


                float vx = 0, 
                      vy = 0,
                      vxOld = 0,
                      vyOld = 0;

                for (int jj = fj; jj < fj + 2; ++jj) {
                    for (int ii = i; ii < i + 2; ++ii) {
                        vx += gridVel[CalcIndexX(ii, jj)] * InterpolateX(ii, jj, particles[k].pos);
                        vxOld += gridVelSaved[CalcIndexX(ii, jj)] * InterpolateX(ii, jj, particles[k].pos);
                    }
                }

                for (int jj = j; jj < j + 2; ++jj) {
                    for (int ii = fi; ii < fi + 2; ++ii) {
                        vy += gridVel[CalcIndexY(ii, jj)] * InterpolateY(ii, jj, particles[k].pos);
                        vyOld += gridVelSaved[CalcIndexY(ii, jj)] * InterpolateY(ii, jj, particles[k].pos);
                    }
                }

                // 速度の更新
                particles[k].vel.x = (1 - flip) * vx + flip * (particles[k].vel.x + vx - vxOld);
                particles[k].vel.y = (1 - flip) * vy + flip * (particles[k].vel.y + vy - vyOld);

                // 壁衝突判定
                // 更新した速度で動かすわけではない！
                particles[k].pos.x = Mathf.Min(Mathf.Max(particles[k].pos.x + vx * dt, 1.001f), n - 1.001f);
                particles[k].pos.y = Mathf.Min(Mathf.Max(particles[k].pos.y + vy * dt, 1.001f), n - 1.001f);
            }
        }

        private void Render() {
            for (int i = 0; i < numParticles; i++) {
                v[i] = new Vector3(particles[i].pos.x, particles[i].pos.y, 0);
                vIdx[i] = i;
            }
            mesh.vertices = v;
            mesh.SetIndices(vIdx, MeshTopology.Points, 0);
            Graphics.DrawMesh(mesh, Matrix4x4.identity, renderMat, 0);
        }

        #region Util
        float InterpolateX(int cellX, int cellY, Vector2 pos) {
            return (1 - Mathf.Abs(cellX - pos.x)) * (1 - Mathf.Abs(cellY + 0.5f - pos.y));
        }

        float InterpolateY(int cellX, int cellY, Vector2 pos) {
            return (1 - Mathf.Abs(cellX + 0.5f - pos.x)) * (1 - Mathf.Abs(cellY - pos.y));
        }

        int CalcIndexX(int i, int j) {

            // |----------|
            // |          |
            // +          +
            // |          |
            // |----------|

            return (i + j * (n + 1));
        }

        int CalcIndexY(int i, int j) {

            // |-----+-----|
            // |           |
            // |           |
            // |           |
            // |-----+-----|

            return y + i + j * n;
        }
        #endregion Util

        #endregion PrivateFuncs
    }

}