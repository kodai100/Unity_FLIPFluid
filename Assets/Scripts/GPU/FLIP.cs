using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;

namespace Kodai.FLIP.GPU {

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
        public ComputeShader FlipCS;

        public Color gizmoColor;
        public Material renderMat;


        private int y;

        [SerializeField] private int numParticles;

        FluidParticle[] particles;
        private float[] gridVel;
        private float[] gridVelSaved;
        private float[] gridMass;
        private float[] gridPressure;
        private float[] gridVelDivergence;
        private int[]   gridFlag;


        static readonly int THREAD_NUM_X = 16;
        ComputeBuffer particlesBuffer;
        ComputeBuffer gridVelocityBuffer;
        ComputeBuffer gridVelocitySavedBuffer;
        ComputeBuffer gridMassBuffer;
        ComputeBuffer gridPressureBuffer;
        ComputeBuffer gridDivergenceBuffer;
        ComputeBuffer gridFlagBuffer;

        // Use this for initialization
        void Start() {
            InitializeGrid();
            InitializeParticle();
            InitializeBuffers();
            SetCSVariables();
        }

        // Update is called once per frame
        void Update() {
            Transfer();
        }

        private void OnDrawGizmos() {

            Gizmos.color = gizmoColor;

            Vector2 a = new Vector2(n, n);
            Gizmos.DrawWireCube(a / 2, a);

            for (int i = 0; i < n; i++) {
                Gizmos.DrawLine(new Vector3(i, 0, 0), new Vector3(i, n, 0));
            }

            for (int i = 0; i < n; i++) {
                Gizmos.DrawLine(new Vector3(0, i, 0), new Vector3(n, i, 0));
            }
        }

        private void OnRenderObject() {
            renderMat.SetPass(0);
            renderMat.SetBuffer("_Particles", particlesBuffer);
            Graphics.DrawProcedural(MeshTopology.Points, numParticles);
        }

        private void OnDestroy() {
            ReleaseBuffer(particlesBuffer);
            ReleaseBuffer(gridVelocityBuffer);
            ReleaseBuffer(gridVelocitySavedBuffer);
            ReleaseBuffer(gridMassBuffer);
            ReleaseBuffer(gridPressureBuffer);
            ReleaseBuffer(gridDivergenceBuffer);
            ReleaseBuffer(gridFlagBuffer);
        }

        void ReleaseBuffer(ComputeBuffer buffer) {
            if(buffer != null) {
                buffer.Release();
            }
        }

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

        void InitializeBuffers() {
            particlesBuffer         = new ComputeBuffer(particles.Length, Marshal.SizeOf(typeof(FluidParticle)));
            gridVelocityBuffer      = new ComputeBuffer(gridVel.Length, Marshal.SizeOf(typeof(float)));
            gridVelocitySavedBuffer = new ComputeBuffer(gridVelSaved.Length, Marshal.SizeOf(typeof(float)));
            gridMassBuffer          = new ComputeBuffer(gridMass.Length, Marshal.SizeOf(typeof(int)));
            gridPressureBuffer      = new ComputeBuffer(gridPressure.Length, Marshal.SizeOf(typeof(float)));
            gridDivergenceBuffer    = new ComputeBuffer(gridVelDivergence.Length, Marshal.SizeOf(typeof(float)));
            gridFlagBuffer          = new ComputeBuffer(gridFlag.Length, Marshal.SizeOf(typeof(int)));

            particlesBuffer.SetData(particles);
        }

        void SetCSVariables() {
            FlipCS.SetInt("_N", n);
            FlipCS.SetFloat("_DT", dt);
            FlipCS.SetFloat("_Flip", flip);
            FlipCS.SetInt("_Y", y);
        }

        void DebugBuffer<T>(int num, ComputeBuffer buffer) {
            T[] bufferData = new T[num];
            buffer.GetData(bufferData);
            
            string str = "GPU: ";
            foreach (T current in bufferData) {
                str += current + ", ";
            }
            Debug.Log(str);
        }

        void DebugParticleBuffer(int num, ComputeBuffer buffer) {
            FluidParticle[] bufferData = new FluidParticle[num];
            buffer.GetData(bufferData);
            string str = "";
            foreach (FluidParticle current in bufferData) {
                str += current.pos + ", ";
            }
            Debug.Log(str);
        }


        /// <summary>
        /// OK
        /// </summary>
        void ResetGrid() {

            int kernel = FlipCS.FindKernel("ResetStaggeredGridQuantities");
            FlipCS.SetBuffer(kernel, "_GridVel", gridVelocityBuffer);
            FlipCS.SetBuffer(kernel, "_GridVelSaved", gridVelocitySavedBuffer);
            FlipCS.SetBuffer(kernel, "_GridMass", gridMassBuffer);
            FlipCS.Dispatch(kernel, gridVel.Length / THREAD_NUM_X, 1, 1);

            kernel = FlipCS.FindKernel("ResetGridQuantities");
            FlipCS.SetBuffer(kernel, "_GridPressure", gridPressureBuffer);
            FlipCS.SetBuffer(kernel, "_GridDivergence", gridDivergenceBuffer);
            FlipCS.SetBuffer(kernel, "_GridFlag", gridFlagBuffer);
            FlipCS.Dispatch(kernel, gridPressure.Length / THREAD_NUM_X, 1, 1);
        }


        void Transfer() {

            ResetGrid();

            int kernel = FlipCS.FindKernel("Transfer");
            FlipCS.SetBuffer(kernel, "_Particles", particlesBuffer);
            FlipCS.SetBuffer(kernel, "_GridVel", gridVelocityBuffer);
            FlipCS.SetBuffer(kernel, "_GridFlag", gridFlagBuffer);
            FlipCS.SetBuffer(kernel, "_GridMass", gridMassBuffer);
            FlipCS.Dispatch(kernel, numParticles / THREAD_NUM_X, 1, 1);
            
            DebugBuffer<float>(gridVel.Length, gridMassBuffer);
            DebugBuffer<int>(gridFlag.Length, gridFlagBuffer);

            //for (int k = 0; k < 2 * n * (n + 1); ++k) {
            //    // 格子に質量があれば
            //    if (gridMass[k] > 1e-8) {
            //        gridMass[k] = 1 / gridMass[k];   // 格子の質量を逆数にして
            //        gridVel[k] *= gridMass[k];       // 速度に掛ける
            //    }
            //}


            //// Save current grid velocities
            //for (int i = 0; i < gridVel.Length; i++) {
            //    gridVelSaved[i] = gridVel[i];
            //}

        }



    }
}