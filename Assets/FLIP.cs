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
        public float flip = 0.96f;
        
        int y;

        private float[] u;  // 半分より手前がx, 後ろがyの値
        private Vector2[] gridVel;  //

        private float[] ux;
        private float[] m;
        private float[] p;
        private float[] divu;

        private Vector2[] pos;
        private Vector2[] particleVel;

        private Cell[] flags;

        private int numParticles;

        private Mesh mesh;
        public Material renderMat;

        void Start() {
            Initialize();

            mesh = new Mesh();
            GetComponent<MeshFilter>().mesh = mesh;
        }
        

        void Update() {
            particles2grid();

            for (int k = 0; k < n * (n + 1); ++k) {
                if (m[y + k] > 1e-8) u[y + k] += -9.8f * dt * (n - 2);
            }

            pressureSolve();
            updateAndAdvect();

            Render();
            
        }

        private void Render() {
            Vector3[] v = new Vector3[numParticles];
            int[] idx = new int[numParticles];
            for (int i = 0; i < numParticles; i++) {
                v[i] = new Vector3(pos[i].x, pos[i].y, 0);
                idx[i] = i;
            }
            mesh.vertices = v;
            mesh.SetIndices(idx, MeshTopology.Points, 0);
            Graphics.DrawMesh(mesh, Matrix4x4.identity, renderMat, 0);
        }

        void Initialize() {
            y = n * (n + 1);    // y軸方向の格子点はx軸方向の格子点より1つ多い
            u = new float[2 * n * (n + 1)];
            ux = new float[2 * n * (n + 1)];
            m = new float[2 * n * (n + 1)];
            p = new float[n * n];
            divu = new float[n * n];
            flags = new Cell[n * n];


            numParticles = (n - 2) / 4 * (n - 4) * 4;
            Debug.Log(numParticles);

            pos = new Vector2[numParticles];
            particleVel = new Vector2[numParticles];

            int idx = 0;
            for (int j = 1; j < n - 1; ++j) {
                for (int i = 1; i < n - 1; ++i) {
                    if (idx < numParticles && i - 1 < (n - 2) / 4) {
                        for (int jj = 0; jj < 2; ++jj) {
                            for (int ii = 0; ii < 2; ++ii) {
                                pos[idx].x = i + 0.25f + ii * 0.5f;
                                pos[idx++].y = j + 0.25f + jj * 0.5f;
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

        void particles2grid() {

            // Initialize Grid Data
            for(int i = 0; i < u.Length; i++) {
                u[i] = 0;
                m[i] = 0;
            }

            // 端の場合SOLID, それ以外の場合EMPTY
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    flags[i + j * n] = (i == 0 || j == 0 || i == n - 1 || j == n - 1) ? Cell.SOLID : Cell.EMPTY;
                }
            }
            
            // 各パーティクルについて
            for (int k = 0; k < numParticles; ++k) {
                // 格子インデックスを取得(格子の大きさは1)
                // TODO 格子の大きさを可変に
                int i  = (int) pos[k].x, 
                    j  = (int) pos[k].y, 
                    fi = (int) (pos[k].x - 0.5f), // 半分左にずらす
                    fj = (int) (pos[k].y - 0.5f); // 半分下にずらす

                flags[i + j * n] = Cell.FLUID;  // 粒子の存在する格子は流体フラグ

                // 上下左右のずれた格子点を参照(x速度)
                for (int jj = fj; jj < fj + 2; ++jj) {
                    for (int ii = i; ii < i + 2; ++ii) {
                        u[IDXX(ii, jj)] += particleVel[k].x * InterpolateX(ii, jj, pos[k]);
                        m[IDXX(ii, jj)] += InterpolateX(ii, jj, pos[k]);
                    }
                }

                // 上下左右の格子点を参照(y速度)
                for (int jj = j; jj < j + 2; ++jj) {
                    for (int ii = fi; ii < fi + 2; ++ii) {
                        u[IDXY(ii, jj)] += particleVel[k].y * InterpolateY(ii, jj, pos[k]);
                        m[IDXY(ii, jj)] += InterpolateY(ii, jj, pos[k]);    // 2回やるの？？
                    }
                }
            }

            for (int k = 0; k < 2 * n * (n + 1); ++k) {
                // 格子に質量があれば
                if (m[k] > 1e-8) {
                    m[k] = 1 / m[k];    // 格子の質量を逆数にして
                    u[k] *= m[k];       // 速度に掛ける
                }
            }

            for (int i = 0; i < u.Length; i++) {
                ux[i] = u[i];
            }
        }

        void pressureSolve() {

            // 格子のi個内側までに関して、速度のダイバージェンスを求める
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

            // Initialize
            for(int i = 0; i < p.Length; i++) {
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

        void updateAndAdvect() {

            // 各粒子について
            for (int k = 0; k < numParticles; ++k) {

                // 自身の格子インデックス
                int i  = (int) (pos[k].x),
                    j  = (int) (pos[k].y), 
                    fi = (int) (pos[k].x - 0.5f), 
                    fj = (int) (pos[k].y - 0.5f);


                float vx = 0, 
                      vy = 0,
                      vxOld = 0,
                      vyOld = 0;

                for (int jj = fj; jj < fj + 2; ++jj) {
                    for (int ii = i; ii < i + 2; ++ii) {
                        vx += u[IDXX(ii, jj)] * (1 - Mathf.Abs(ii - pos[k].x)) * (1 - Mathf.Abs(jj + 0.5f - pos[k].y));
                        vxOld += ux[IDXX(ii, jj)] * (1 - Mathf.Abs(ii - pos[k].x)) * (1 - Mathf.Abs(jj + 0.5f - pos[k].y));
                    }
                }

                for (int jj = j; jj < j + 2; ++jj) {
                    for (int ii = fi; ii < fi + 2; ++ii) {
                        vy += u[IDXY(ii, jj)] * (1 - Mathf.Abs(ii + 0.5f - pos[k].x)) * (1 - Mathf.Abs(jj - pos[k].y));
                        vyOld += ux[IDXY(ii, jj)] * (1 - Mathf.Abs(ii + 0.5f - pos[k].x)) * (1 - Mathf.Abs(jj - pos[k].y));
                    }
                }

                // 速度の更新
                particleVel[k].x = (1 - flip) * vx + flip * (particleVel[k].x + vx - vxOld);
                particleVel[k].y = (1 - flip) * vy + flip * (particleVel[k].y + vy - vyOld);

                // 壁衝突判定
                // 更新した速度で動かすわけではない！
                pos[k].x = Mathf.Min(Mathf.Max(pos[k].x + vx * dt, 1.001f), n - 1.001f);
                pos[k].y = Mathf.Min(Mathf.Max(pos[k].y + vy * dt, 1.001f), n - 1.001f);
            }
        }

        int IDXX(int i, int j) {
            return (i + j * (n + 1));
        }

        int IDXY(int i, int j) {
            return y + i + j * n;
        }

    }

}