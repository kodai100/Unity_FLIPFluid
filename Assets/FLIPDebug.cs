using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Kodai.FLIP.CPU {

    public class FLIPDebug : MonoBehaviour {

        public enum Cell {
            EMPTY, FLUID, SOLID
        }

        struct FluidParticle {
            public Vector2 pos;
            public Vector2 vel;
        }
        FluidParticle[] particles;

        public int n = 32;
        public float dt = 1 / 120f;
        public float flip = 0.96f;
        

        private float[] gridVelX;
        private float[] gridVelY;

        private float[] ux;
        private float[] uy;

        private float[] mx;
        private float[] my;

        private float[] p;
        private float[] divu;
        private Cell[] flags;

        private int numParticles;

        private Mesh mesh;
        public Material renderMat;

        Texture2D texU;
        Texture2D texDiv;
        Texture2D texFlag;

        void Start() {
            Initialize();

            mesh = new Mesh();
            GetComponent<MeshFilter>().mesh = mesh;

            texU = new Texture2D((n+1) * n * 2, 1);
            texU.filterMode = FilterMode.Point;

            texDiv = new Texture2D(n * n, 1);
            texDiv.filterMode = FilterMode.Point;

            texFlag = new Texture2D(n * n, 1);
            texFlag.filterMode = FilterMode.Point;
        }


        void Update() {
            particles2grid();

            // y軸方向に重力を作用させる
            for (int k = 0; k < n * (n + 1); ++k) {
                if (my[k] > 1e-8) gridVelY[k] += -9.8f * dt * (n - 2);  // (n - 2)はシミュレーションスケール？
            }

            pressureSolve();
            updateAndAdvect();

            Render();
            TextureRender();
        }

        void VelTex() {
            for (int i = 0; i < gridVelX.Length + uy.Length; i++) {

                if (gridVelX.Length > i) {
                    texU.SetPixel(i, 0, new Color(gridVelX[i], gridVelX[i], gridVelX[i]));
                } else {
                    texU.SetPixel(i, 0, new Color(gridVelY[i - gridVelX.Length], gridVelY[i - gridVelX.Length], gridVelY[i - gridVelX.Length]));
                }

            }
            texU.Apply();
        }

        void TextureRender() {
            

            for (int i = 0; i < divu.Length; i++) {
                texDiv.SetPixel(i, 0, new Color(divu[i], divu[i], divu[i]));

                if (flags[i] == Cell.FLUID) {
                    texFlag.SetPixel(i, 0, Color.blue);
                } else if (flags[i] == Cell.SOLID) {
                    texFlag.SetPixel(i, 0, Color.green);
                } else {
                    texFlag.SetPixel(i, 0, Color.black);
                }
            }
            texDiv.Apply();
            texFlag.Apply();

        }

        private void OnGUI() {
            GUI.DrawTexture(new Rect(0, 10, Screen.width, 10), texDiv);
            GUI.DrawTexture(new Rect(0, 50, Screen.width, 10), texFlag);

            GUI.DrawTexture(new Rect(0, 80, Screen.width, 10), texU);
        }

        private void OnDrawGizmos() {
            Vector2 a = new Vector2(n, n);
            Gizmos.DrawWireCube(a / 2, a);
        }

        private void Render() {
            Vector3[] v = new Vector3[numParticles];
            int[] idx = new int[numParticles];
            for (int i = 0; i < numParticles; i++) {
                v[i] = new Vector3(particles[i].pos.x, particles[i].pos.y, 0);
                idx[i] = i;
            }
            mesh.vertices = v;
            mesh.SetIndices(idx, MeshTopology.Points, 0);
            Graphics.DrawMesh(mesh, Matrix4x4.identity, renderMat, 0);
        }

        void Initialize() {
            gridVelX = new float[n * (n + 1)];
            gridVelY = new float[n * (n + 1)];
            ux = new float[n * (n + 1)];    // v
            uy = new float[n * (n + 1)];
            mx = new float[n * (n + 1)];
            my = new float[n * (n + 1)];

            p = new float[n * n];
            divu = new float[n * n];
            flags = new Cell[n * n];


            numParticles = (n - 2) / 4 * (n - 4) * 4;
            particles = new FluidParticle[numParticles];
            Debug.Log(numParticles);

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

        float InterpolateX(int cellX, int cellY, Vector2 pos) {
            return (1 - Mathf.Abs(cellX - pos.x)) * (1 - Mathf.Abs(cellY + 0.5f - pos.y));
        }

        float InterpolateY(int cellX, int cellY, Vector2 pos) {
            return (1 - Mathf.Abs(cellX + 0.5f - pos.x)) * (1 - Mathf.Abs(cellY - pos.y));
        }

        void InitializeGrid() {
            // Initialize Grid Data
            for (int i = 0; i < gridVelX.Length; i++) {
                gridVelX[i] = 0;
                gridVelY[i] = 0;
                mx[i] = 0;
                my[i] = 0;
            }


            // 端の場合SOLID, それ以外の場合EMPTY
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    flags[i + j * n] = (i == 0 || j == 0 || i == n - 1 || j == n - 1) ? Cell.SOLID : Cell.EMPTY;
                }
            }
        }

        void particles2grid() {

            InitializeGrid();

            // 各パーティクルについて
            for (int k = 0; k < numParticles; ++k) {
                // パーティクルの属する格子インデックスを取得(格子の大きさは1)
                // TODO 格子の大きさを可変に
                int i = (int)particles[k].pos.x,
                    j = (int)particles[k].pos.y,
                    fi = (int)(particles[k].pos.x - 0.5f), // 半分左にずらす
                    fj = (int)(particles[k].pos.y - 0.5f); // 半分下にずらす

                flags[i + j * n] = Cell.FLUID;  // 粒子の存在する格子は流体フラグ


                // 上下左右のずれた格子点を参照(x速度)
                for (int jj = fj; jj < fj + 2; ++jj) {
                    for (int ii = i; ii < i + 2; ++ii) {
                        gridVelX[IDXX(ii, jj)] += particles[k].vel.x * InterpolateX(ii, jj, particles[k].pos);
                        mx[IDXX(ii, jj)] += InterpolateX(ii, jj, particles[k].pos);
                    }
                }

                // 上下左右の格子点を参照(y速度)
                for (int jj = j; jj < j + 2; ++jj) {
                    for (int ii = fi; ii < fi + 2; ++ii) {
                        
                        gridVelY[IDXY(ii, jj)] += particles[k].vel.y * InterpolateY(ii, jj, particles[k].pos);
                        my[IDXY(ii, jj)] += InterpolateY(ii, jj, particles[k].pos);
                    }
                }
            }
            
            for (int k = 0; k < n * (n + 1); ++k) {
                // 格子に質量があれば
                if (mx[k] > 1e-8f) {
                    mx[k] = 1 / mx[k];    // 格子の質量を逆数にして
                    gridVelX[k] *= mx[k];       // 速度に掛ける
                }

                if(my[k] > 1e-8f) {
                    my[k] = 1 / my[k];    // 格子の質量を逆数にして
                    gridVelY[k] *= my[k];       // 速度に掛ける
                }
            }

            VelTex();

            for (int i = 0; i < gridVelX.Length; i++) {
                ux[i] = gridVelX[i];
                uy[i] = gridVelY[i];
            }
        }

        void pressureSolve() {
            
            // 格子のi個内側までに関して、速度のダイバージェンスを求める
            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    divu[i + j * n] = -(gridVelX[i + 1 + j * (n + 1)] - gridVelX[i + j * (n + 1)] + gridVelY[i + (j + 1) * n] - gridVelY[i + j * n]);
                }
            }

            // 境界条件
            for (int k = 1; k < n - 1; ++k) {
                divu[1 + k * n] -= gridVelX[1 + k * (n + 1)];   // 左端
                divu[n - 2 + k * n] += gridVelX[n - 1 + k * (n + 1)]; // 右端
                divu[k + n] -= gridVelY[k + n];
                divu[k + (n - 2) * n] += gridVelY[k + (n - 1) * n];
            }

            // Initialize
            for (int i = 0; i < p.Length; i++) {
                p[i] = 0;
            }

            // ガウスザイデル法でラプラス方程式を解く
            for (int k = 0; k < 10 * n * n; ++k) {    // Gauss-Seidel solve for pressure
                // 格子の1つ内側について
                for (int j = 1; j < n - 1; ++j) {
                    for (int i = 1; i < n - 1; ++i) {

                        // 流体以外であれば処理をスルー
                        if ((int)flags[i + j * n] != 1) continue;

                        // 上下左右を判定し、枠だった場合は0
                        float l = i > 1 ? -1 : 0,   // left
                              r = i < n - 2 ? -1 : 0,   // right
                              b = j > 1 ? -1 : 0,       // bottom
                              t = j < n - 2 ? -1 : 0,   // top
                              c = -(l + r + b + t);     // center

                        // 圧力の計算
                        p[i + j * n] = 1 / c * (divu[i + j * n] - l * p[i - 1 + j * n] - r * p[i + 1 + j * n] - b * p[i + j * n - n] - t * p[i + j * n + n]);
                    }
                }
            }


            // 圧力による速度の更新

            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n; ++i) {
                    gridVelX[i + j * (n + 1)] -= p[i + j * n] - p[i - 1 + j * n];
                }
            }

            for (int j = 1; j < n; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    gridVelY[i + j * n] -= p[i + j * n] - p[i + j * n - n];
                }
            }

            // 端を0に
            for (int k = 1; k < n - 1; ++k) {
                gridVelX[1 + k * (n + 1)] = 0;
                gridVelX[n - 1 + k * (n + 1)] = 0;
                gridVelY[k + n] = 0;
                gridVelY[k + (n - 1) * n] = 0;
            }
        }

        void updateAndAdvect() {

            // 各粒子について
            for (int k = 0; k < numParticles; ++k) {

                // 自身の格子インデックス
                int i = (int)(particles[k].pos.x),
                    j = (int)(particles[k].pos.y),
                    fi = (int)(particles[k].pos.x - 0.5f),
                    fj = (int)(particles[k].pos.y - 0.5f);


                float vx = 0,
                      vy = 0,
                      vxOld = 0,
                      vyOld = 0;

                for (int jj = fj; jj < fj + 2; ++jj) {
                    for (int ii = i; ii < i + 2; ++ii) {
                        vx += gridVelX[IDXX(ii, jj)] * InterpolateX(ii, jj, particles[k].pos);
                        vxOld += ux[IDXX(ii, jj)] * InterpolateX(ii, jj, particles[k].pos);
                    }
                }
                if (k == 0) Debug.Log(vx);

                for (int jj = j; jj < j + 2; ++jj) {
                    for (int ii = fi; ii < fi + 2; ++ii) {
                        vy += gridVelY[IDXY(ii, jj)] * InterpolateY(ii, jj, particles[k].pos);
                        vyOld += ux[IDXY(ii, jj)] * InterpolateY(ii, jj, particles[k].pos);
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

        int IDXX(int i, int j) {
            return (i + j * (n + 1));
        }

        int IDXY(int i, int j) {
            return i + j * n;
        }

    }

}