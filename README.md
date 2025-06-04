**_To fit discrete extended generalized Pareto distribution (degpd) and zero-inflated discrete extended generalized Pareto distribution (zidegpd), we developed a code for new families and run it using evgam package fucntion_**

**_Intallation_**
1. Download the code and set the working directory. e.g., setwd("C:\Users\atouqeer\Downloads\degpd-and-zidegpd-main\degpd-and-zidegpd-main")
2. Call all the C++ and R functions using the following code

```markdown
dest <- "./R/"      # this function all R function 
files = list.files(dest, full.names = T)
for (i in 1:length(files)) {
  source(files[i])
}


dest <- "./src/"  
files = list.files(dest, full.names = T)
for (i in 1:length(files)) {
  Rcpp::sourceCpp(files[i])
}
```





### Fitting of degpd-and-zidegpd using simulation example

We run an example on the simulation data to show how the code works. The "Fit_degpd_zidegpd.R" file simulated the data and fit simply the discrete extended generalized Pareto distribution (degpd) and zero-inflated discrete extended generalized Pareto distribution (zidegpd) with their GAM forms as well. These models are proposed in the paper "Ahmad T, Gaetan C, Naveau P. An extended generalized Pareto regression model for count data. Statistical Modelling. 2024;0(0). [doi:10.1177/1471082X241266729](https://journals.sagepub.com/doi/abs/10.1177/1471082X241266729)".

We are using the functions of evgam R package (Youngman, 2020): An R package for Generalized Additive Extreme Value Models. 
https://doi.org/10.48550/arXiv.2003.04067 behind to run our own developed R code.

The example with fitting of the degpd 1 model is shown in the "Fit_degpd_zidegpd.R" code. The other degpd models, 2, 3, and 4, can be fitted by changing **m**. 

---
In the code, the m=1 is corresponding to model $$G\left(u; \psi\right)={u}^{\kappa},$$

---
m=2 corresponds to the model
$$G\left(u;\psi\right)= p{u}^{\kappa_1} + \left(1-p\right){u}^{\kappa_2},$$

---
m=3 corresponds to the model
$$G\left(u;\psi\right)=1-D_{\delta}\{\left(1-u\right)^{\delta}\},$$

---

and m=4 corresponds to the model
$$G\left(u;\psi\right)=\left[1-D_{\delta}\{(1-u)^{\delta}\}\right]^{\kappa/2}$$

---
The zidegpd models can also be fitted by changing family " **degpd**" to "**zidegpd**" and by changing **m**. 

The m=1 is corresponding to model $$G\left(u; \psi\right)={u}^{\kappa},$$


$$m=2$$ is corresponding to model (**not developed yet**)
$$G\left(u;\psi\right)= p{u}^{\kappa_1} + \left(1-p\right){u}^{\kappa_2},$$


m=3 corresponds to the model
$$G\left(u;\psi\right)=1-D_{\delta}\{\left(1-u\right)^{\delta}\},$$


and m=4 corresponds to the model
$$G\left(u;\psi\right)=\left[1-D_{\delta}\{(1-u)^{\delta}\}\right]^{\kappa/2}$$


**Note** 
1. Ignore the error
"Error in Rcpp::sourceCpp(files[i]) : 
  The filename 'Makevars' does not have an extension of .cc or .cpp, so it cannot be compiled." The code is working correctly; this error appears from the evgam package code. We are using some functions behind the scenes from the evgam package. We will surely introduce our proposed models as new families in evgam.
  
2. We developed and executed this code on the Windows operating system. You may face problems when you run on a macOS (some instructions may be helpful: https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/). OR, follow the instructions below.


# If you are using macOS (especially on Apple Silicon) and want to compile C++ files with `Rcpp`, follow these steps to set up `gfortran` and avoid linker errors.

## 🛠 Setup: GFortran + Rcpp on macOS (Apple Silicon)

---

### Step 1: Reinstall GCC via Homebrew

```bash
brew reinstall gcc

```
### Step 2: Verify Installed GFortran Version
```bash
ls /opt/homebrew/opt/gcc/bin | grep gfortran

```
Check the version (replace with your installed version if it differs): In my case, I have version 15 as gfortran-15 --version

---
### Step 3: Configure R to Use GFortran
Edit (or create) the file ~/.R/Makevars and add:
```bash
nano ~/.R/Makevars
```
Paste the following (replace 15 with your version if different):
```
FC = /opt/homebrew/opt/gcc/bin/gfortran-15
F77 = /opt/homebrew/opt/gcc/bin/gfortran-15

FLIBS = -L/opt/homebrew/opt/gcc/lib/gcc/15 -lgfortran -lquadmath -lm

CXXFLAGS = -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/opt/gcc/lib/gcc/15 -lgfortran -lquadmath -lm

```
Save (Ctrl+O) the file and exit (Ctrl+X) the editor.

---
### Step 4: Update Your Terminal PATH
Make sure your terminal knows where to find gfortran:
```
echo 'export PATH="/opt/homebrew/opt/gcc/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

---
### Step 5: Step 5: Compile Your C++ Files in R
Open R and run:

```
system("gfortran-15 --version")

dest <- "./src/"
files <- list.files(dest, full.names = TRUE)
for (i in seq_along(files)) {
  Rcpp::sourceCpp(files[i])
}
```

---
**For further discussion, feel free to write me at (touqeer.ahmad8960@gmail.com)**
 

