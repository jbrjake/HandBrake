﻿/*  AudioEncoder.cs $
    This file is part of the HandBrake source code.
    Homepage: <http://handbrake.fr>.
    It may be used under the terms of the GNU General Public License. */

namespace HandBrake.ApplicationServices.Model.Encoding
{
    using System.ComponentModel;

    public enum AudioEncoder
    {
        [Description("AAC (faac)")]
        Faac = 0,

        [Description("MP3 (lame)")]
        Lame,

        [Description("AC3 (ffmpeg)")]
        Ac3,

        [Description("Passthrough (AC3/DTS)")]
        Passthrough,

        [Description("Passthrough (AC3)")]
        Ac3Passthrough,

        [Description("Passthrough (DTS)")]
        DtsPassthrough,

        [Description("Vorbis (vorbis)")]
        Vorbis
    }
}