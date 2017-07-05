using System;
using UnityEngine;


namespace UTJ.FrameCapturer
{
    public class WebMContext : MovieRecorderContext
    {
        [Serializable]
        public class EncoderConfig
        {
            [HideInInspector] public bool captureVideo = true;
            [HideInInspector] public bool captureAudio = true;
            public int videoBitrate = 8192000;
            public int audioBitrate = 64000;
        }

        fcAPI.fcWebMContext m_ctx;
        fcAPI.fcStream m_ostream;
        EncoderConfig m_config;

        public override Type type { get { return Type.WebM; } }

        public override void Initialize(MovieRecorder recorder)
        {
            m_config = recorder.webmConfig;
            fcAPI.fcWebMConfig webmconf = fcAPI.fcWebMConfig.default_value;
            webmconf = fcAPI.fcWebMConfig.default_value;
            webmconf.video = m_config.captureVideo;
            webmconf.audio = m_config.captureAudio;
            webmconf.video_width = recorder.scratchBuffer.width;
            webmconf.video_height = recorder.scratchBuffer.height;
            webmconf.video_target_framerate = 60;
            webmconf.video_target_bitrate = m_config.videoBitrate;
            webmconf.audio_target_bitrate = m_config.audioBitrate;
            webmconf.audio_sample_rate = AudioSettings.outputSampleRate;
            webmconf.audio_num_channels = fcAPI.fcGetNumAudioChannels();
            m_ctx = fcAPI.fcWebMCreateContext(ref webmconf);

            var path = recorder.outputDir.GetFullPath() + "/" + DateTime.Now.ToString("yyyyMMdd_HHmmss") + ".webm";
            m_ostream = fcAPI.fcCreateFileStream(path);
            fcAPI.fcWebMAddOutputStream(m_ctx, m_ostream);
        }

        public override void Release()
        {
            fcAPI.fcGuard(() =>
            {
                m_ctx.Release();
                m_ostream.Release();
            });
        }

        public override void AddVideoFrame(byte[] frame, fcAPI.fcPixelFormat format, double timestamp)
        {
            if (m_config.captureVideo)
            {
                fcAPI.fcWebMAddVideoFramePixels(m_ctx, frame, format, timestamp);
            }
        }

        public override void AddAudioFrame(float[] samples, double timestamp)
        {
            if (m_config.captureAudio)
            {
                fcAPI.fcWebMAddAudioFrame(m_ctx, samples, samples.Length);
            }
        }
    }
}
