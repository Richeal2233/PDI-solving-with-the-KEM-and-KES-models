# PDI-solving-with-the-KEM-and-KES-models
A program of calculating the growth rate of parametric instabilities (PDI) of RF waves with the kinetic electromagnetic model (KEM) and the kinetic electrostatic model(KES)

你好！欢迎使用我们基于matlab计算参数不稳定性（PDI）的全动理学电磁（KEM）模型和全动理学静电（KES）模型程序！具体如下。

求解过程的核心文件是solvew.m。第一部分是所有等离子体参数和波参数的输入和归一化处理。
第二部分是matlab上使用patternsearch算法需要设置的初始值“wiv”和求解选项“PSoptions”、lb（下限）、ub（上限）等，需要此算法在matlab中已安装。
第三部分是求解目标函数fw(x)并输出解（由1e9归一化的低频波的实际频率和增长率）。

目标函数为fw.m。首先，它接收输入值并将其转换为归一化频率。
并且它根据PIM.m和频率计算色散矩阵（6x6）的条件数的相反数。
在满足矩阵行列式接近零的解条件时，其条件数达到最大值，其相反数（目标函数）达到最小值，这是patternsearch算法的目标。

KEM和KES模型的色散矩阵由PIM.m处理，其中包含线性、准线性和非线性阶次的所有介电张量矩阵。
因此，PIM.m是物理模型和所有剩余文件的核心。它们都是色散矩阵的一部分。您可以在PIM.m中更改模型（KEM或KES）。

整个程序可以通过改变波参数（主要是“nper”，低频波的垂直折射率）来计算衰减通道内某点的增长率。
等离子体参数也可以根据不同的场景，甚至不同的装置进行更改。
但值得注意的是，模型的色散矩阵非常复杂，因此目标函数具有很强的奇异性，计算非常耗时。
此外，初始值的一点变化都可能会导致无法找到解，更不用说波参数和等离子体参数的变化了。
因此，这些模型适用于在已知衰变通道中与一些简化模型进行比对，而不适合在广泛的参数空间中寻找衰变通道。

有一个在LHCD期间PDI的离子回旋准模（ICQM）的示例点，具有JET上SOL的典型参数。
波参数“nper”和初始值“wiv”已经设置好，可以直接运行solvew.m获得结果。请等待几十分钟甚至一个小时，在某些积分中忽略警告。它们对结果并不重要。

程序逻辑可以进一步优化。版权归清华大学及通讯作者所有gaozhe@tsinghua.edu.cn, huangzk23@mails.tsinghua.edu.cn，相关文献https://doi.org/10.1063/1.5139281.


Hello! Welcome to our program of the kinetic electromagnetic (KEM) model and the kinetic electrostatic (KES) model of parametric instabilities (PDI) based on matlab! It's organized by the following method.

The core file of the solving process is solvew.m. 
The first part is the input and normalization treatment of all the plasma and wave parameters.
The second part is the settings of the initial value "wiv" and solving options "PSoptions", lb(lower bound), ub(upper bound) of the patternsearch on matlab, which also need to be installed in matlab.
The third part is to solve the objective function fw(x) and output the solution (the real frequency and the growth rate of the low frequency wave normalized by 1e9).

The objective function is fw.m.
First, it receives the input value and turn it into the normalized frequency.
And it calculate the opposite of the condition number of the dispersion matrix (6x6) based on the PIM.m and the frequency.
At the solution where the matrix determinant close to zero, the condition number reaches a maximum value and its opposite number (objective function) reachs a minimum value, where is the goal of the patternsearch algorithm.

The dispersion matrix of the KEM and KES model is handled by PIM.m, which containing all the dielectric tensor matrix of linear, quasi-linear and nonlinear orders.
So the PIM.m is the core of the physical model and all remaining files. They are all part of the dispersion matrix. You can change the model (KEM or KES) in PIM.m

The whole program can calculate the points of the decay channels by changing the wave parameters (mainly "nper", the perpendicular refractive index of the low frequency wave).
The plasma prameters can also be changed for different scenarios, even different devices.
But it is noted that the dispersion matrix of the models are very complex, so the objective function has strong singularity and the calculation is very time-consuming.
Additionally, a little change of the initial value could lead to the inability to find a solution, let alone changes of the wave parameters and plasma parameters.
So, the models are suitable for the benchmark with some simplified models in the known decay channels,  rather than searching for the decay channels in a broad parameter space.

There is an example point of ion cyclotron quasi mode (ICQM) of PDI during LHCD with the typical parameters of SOL on JET.
The wave parameter "nper", initial value "wiv" are set up and solvew.m can be run directly to get the corresponding results.
Please waiting for several tens of minutes or even an hour and ignoring the warinings during some of the integrals. They are not important for the results.

The program logic can be optimized further. The copyright belongs to Tsinghua University and the corresponding author gaozhe@tsinghua.edu.cn, huangzk23@mails.tsinghua.edu.cn, related literature https://doi.org/10.1063/1.5139281.
