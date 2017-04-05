
<table>
    <thead>
        <tr>
            <th>ID</th>
            <th>Name</th>
            <th>Description</th>
            <th>Comments</th>
        </tr>
    </thead>
    <tbody>
        <tr><td>0</td><td>flat tpdf</td><td>no noise shaping</td><td>no error correction feedback</td></tr>
        <tr><td>1</td><td>classic</td><td>subtle noise shaping (HF emphasis)</td><td>original noise shaping from older versions of ReSampler</td></tr>
        <tr><td>2</td><td>flat tpdf (with error-correction feedback)</td><td>gentle highpass response</td><td>error correction feedback loop results in first-order (6dB/oct) slope</td></tr>
        <tr><td>3</td><td>Modified E-Weighted</td><td>notch in 3-4kHz range (which our ears are sensitive to)</td><td>from the paper "Minimally Audible Noise Shaping"</td></tr>
        <tr><td>4</td><td>Wannamaker 3-tap</td><td>simple f-weighted curve, with notch around 4kHz</td><td>from the paper "Psychoacoustically Optimal Noise Shaping"</td></tr>
        <tr><td>5</td><td>Lipshitz</td><td>e-weighted curve with notches around 4k and 12k</td><td>from the paper "Minimally Audible Noise Shaping"</td></tr>
        <tr><td>6</td><td>standard</td><td>Smooth curve with notches at 3150 and 11250Hz, and extreme HF emphasis</td><td>default noise shape</td></tr>
        <tr><td>7</td><td>Wannamaker 24-tap</td><td>notches around 3.5kHz and 12kHz</td><td>from the paper "Psychoacoustically Optimal Noise Shaping"</td></tr>
        <tr><td>8</td><td>Wannamaker 9-tap</td><td>notches around 3.5kHz and 12kHz</td><td>from the paper "Psychoacoustically Optimal Noise Shaping"</td></tr>
        <tr><td>9</td><td>High28</td><td>notches at 3150Hz and 11.25kHz with 28dB high shelf</td><td></td></tr>
        <tr><td>10</td><td>Improved E-Weighted</td><td>widely used in many DAWs and audio software</td><td>from the paper "Minimally Audible Noise Shaping"</td></tr>
        <tr><td>11</td><td>High30</td><td>notches at 3150Hz and 11.25kHz with 30dB high shelf</td><td>sounds great to older listeners, but high-frequencies may annoy younger listeners. Nevertheless, at the intended playback volume, the high-frequency components will still be very quiet relative to the program material</td></tr>
        <tr><td>12</td><td>High32</td><td>notches at 3150Hz and 11.25kHz with 32dB high shelf</td><td>sounds great to older listeners, but high-frequencies may annoy younger listeners. Nevertheless, at the intended playback volume, the high-frequency components will still be very quiet relative to the program material</td></tr>
    </tbody>
</table>



