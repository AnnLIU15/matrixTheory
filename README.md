# matrix Theory

2022-autumn sysu-seit matrix theory final project

|                          Author                          |                            Email                            |
| :-------------------------------------------------------: | :----------------------------------------------------------: |
| [AnnLIU15 (ZhaoY) (github.com)](https://github.com/AnnLIU15) | [liuzhy86@mail2.sysu.edu.cn](mailto:liuzhy86@mail2.sysu.edu.cn) |
|      [zxh1128 (github.com)](https://github.com/zxh1128)      |                                                              |

## git 教程

* 每次打开项目的时候，先在terminal输入

  ```
  git pull origin master
  ```
* 每次写完一点东西的时候，在terminal输入

  ```
  # 将所有更改提交缓存区
  git add .  
  # 将所有缓存区内的变动提交
  git commit -m "本次commit的名字"
  # 提交到远程
  git push origin master
  ```

## 项目结构

> * algo        -- 核心算法文件
> * utils        -- 除了核心算法外的function（yaml基于[Google Code Archive - Long-term storage for Google Code Project Hosting.](https://code.google.com/archive/p/yamlmatlab/)）
> * assets     -- 图像资源
> * output    -- 算法输出结果

## 选题与需要实现的算法

* Y. Hu, D. Zhang, J. Ye, X. Li, and X. He, “Fast and Accurate Matrix Completion via Truncated Nuclear Norm Regularization,” *Ieee T Pattern Anal*, vol. 35, no. 9, pp. 2117–2130, Sep. 2013.
* https://ieeexplore.ieee.org/document/6389682/?arnumber=6389682
* - [X] TNNR-ADMM
* - [X] TNNR-APGL
* - [ ] TNNR-ADMMAP -- 字符迭代

---

对比算法

- [X] opt-space
- [X] SVT & SVP[HauLiang/Matrix-Completion-Methods: Several Classic Low-Rank Matrix Completion Methods, such as SVP, SVT, Sp-lp, and TNNR-ADMM... (github.com)](https://github.com/HauLiang/Matrix-Completion-Methods)

## 要求

* Write a detailed report in LaTex on your whole project, including the background, problem formulation, mathematical method(s) and solution(s), algorithm(s), simulation setting and results, discussions, and conclusions;
* Compose a compact PPT on your report and make a presentation to your classmates if asked. The length of the PPT should be no longer than 20 slides;
* Submit your report, PPT, and the source codes of your simulation experiments to the TA no later than December 9th, 2022.
