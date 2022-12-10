# matrix Theory

2022-autumn sysu-seit matrix theory final project

|                          Author                          |                             Email                             |
| :-------------------------------------------------------: | :------------------------------------------------------------: |
| [AnnLIU15 (ZhaoY) (github.com)](https://github.com/AnnLIU15) |  [liuzhy86@mail2.sysu.edu.cn](mailto:liuzhy86@mail2.sysu.edu.cn)  |
|      [zxh1128 (github.com)](https://github.com/zxh1128)      | [zhengxh56@mail2.sysu.edu.cn](mailto:zhengxh56@mail2.sysu.edu.cn) |

<div align='center'>
    <font color="red" size=10><b>华哥 我的超人 &#x1F645&#x200D&#x2640&#xFE0F;</b></font>
</div>


<div align='center'>
    <font color="red" size=10><b>炀神 真正的超人&#x1F9B8;</b></font>
</div>
~~点解html无法在GitHub上正常显示字号和颜色捏~~

## 报告与PPT

* 报告与PPT的 `svg`与 `eps`图片都使用[Inkscpae](www.inkscape.org/)软件进行**边缘裁剪**
* 报告源文件 [./assets/report.tex](./assets/report.tex)
  * 使用的是neurips2022模板
  * pdf为[./assets/report.pdf](./assets/report.pdf)
    * 使用的是xelatex+pdftex
    * vscode+texlive编译，配置文件 [.vscode/settings.json](.vscode/settings.json)需要配合 `LaTeX`与 `LaTeX Workshop`插件使用，编译流程相同的情况下应该可以在**texstudio**等编辑器上运行
* PPT [./assets/tnnr.pptx](./assets/tnnr.pptx)

## git 教程

* 每次打开项目的时候，先在terminal输入

  ```bash
  git pull origin master
  ```
* 每次写完一点东西的时候，在terminal输入

  ```bash
  # 将所有更改提交缓存区
  git add .  
  # 将所有缓存区内的变动提交
  git commit -m "本次commit的名字"
  # 提交到远程
  git push origin master
  ```

## 运行方式

以fig1和fig9为例子 (需要快速运行请在yaml文件设置为 `min_R = max_R`)

### Fig1

```matlab
addpath(genpath(cd)); run("scripts/getFig1.m")
```

### Fig9

```matlab
addpath(genpath(cd)); run("scripts/getFig9.m")
```

## 项目结构

> * scripts	   -- 运行脚本（由于本人一开始使用vscode运行，故在matlab上运行需要先在下方运行
>
>   ```matlab
>   addpath(genpath(cd));
>   ```
> * algo         -- 核心算法文件
> * utils        -- 除了核心算法外的function
>
>   * yaml基于[Google Code Archive - Long-term storage for Google Code Project Hosting.](https://code.google.com/archive/p/yamlmatlab/)
>   * tight_subplot基于[tight_subplot(Nh, Nw, gap, marg_h, marg_w) - File Exchange - MATLAB Central (mathworks.cn)](https://ww2.mathworks.cn/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
> * assets    -- 图像资源
> * output    -- 算法输出结果
> * report    -- 大作业报告、ppt（tex源文件、pdf以及pptx）

## 选题与需要实现的算法

Baseline/对比算法

- [X] SVT [J.-F. Cai, E. J. Candes, and Z. Shen, “A Singular Value Thresholding Algorithm for Matrix Completion.” arXiv, Oct. 17, 2008. ](http://arxiv.org/abs/0810.3286)
- [X] SVP [P. Jain, R. Meka, and I. Dhillon, “Guaranteed Rank Minimization via Singular Value Projection,” in *Advances in Neural Information Processing Systems*, 2010, vol. 23.](https://proceedings.neurips.cc/paper/2010/hash/08d98638c6fcd194a4b1e6992063e944-Abstract.html)
- [X] OptSpace [R. H. Keshavan, A. Montanari and S. Oh, &#34;Matrix Completion From a Few Entries,&#34; in IEEE Transactions on Information Theory, vol. 56, no. 6, pp. 2980-2998, June 2010](https://ieeexplore.ieee.org/document/5466511/).

---

[Y. Hu, D. Zhang, J. Ye, X. Li, and X. He, “Fast and Accurate Matrix Completion via Truncated Nuclear Norm Regularization,” *Ieee T Pattern Anal*, vol. 35, no. 9, pp. 2117–2130, Sep. 2013.](https://ieeexplore.ieee.org/document/6389682/?arnumber=6389682)

* - [X] TNNR-ADMM
* - [X] TNNR-APGL
* - [X] TNNR-ADMMAP

>  在我们复现的代码中，归一化会有更好的效果，Fig4~8未作归一化（结果更接近），Fig9归一化了（不想调超参了）。复现结果与作者论文展示的效果存在差异，作者参数在本复现代码中无法发挥良好，不清楚作者图像是否做了预处理，以下为一些差异：
>
>  * 运算时间：我们的复原在运算时间上无法完成tnnr-admmap快于tnnr-admm的效果，可能是由于我们设置的$\rho_0$较为保守（按照作者的设置会发散）
>  * 效果：APGL在图像复原的效果无法达到作者展示的清晰度
>  * Baseline 在人工生成环境表现过好，在图像mask环境中表现较差
>
>  由于没时间了，没法调最好的超参Orz

## 要求

* Write a detailed report in LaTex on your whole project, including the background, problem formulation, mathematical method(s) and solution(s), algorithm(s), simulation setting and results, discussions, and conclusions;
* Compose a compact PPT on your report and make a presentation to your classmates if asked. The length of the PPT should be **no longer than 20 slides**;
* Submit your **report, PPT, and the source codes of your simulation experiments** to the TA no later than **December 9th, 2022**.
