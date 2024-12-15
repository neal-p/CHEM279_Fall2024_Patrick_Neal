# REPO LINK: [https://github.com/neal-p/CHEM279](https://github.com/neal-p/CHEM279)

# HW5

## Compilation Instructions
  1. Go to the top level repo directory `CHEM279/`
  2. `mkdir -p build`
  3. `cd build`
  4. `cmake ..`
  5. `make hw5`

## Run Instructions
Once you have compiled, there will be one executable:
  - `CHEM279/HW5/hw5`

There is a script `run_samples.sh` that will build and run the code on all `sample_input` files putting output in a corresponding `my_outputs` folder.
NOTE! I was only successful for `H2` in this problem set :( 

### Main hw5 program
This executable takes in an xyz file in the specification of our class and a basis set file in Gaussian format. It gets the initial energy and computes the gradient. An exmaple output is shown below for `HF`:
`./hw5 sample_input/H2.txt basis/basis_set.txt`

```
Nuclear Repulsion Energy=19.4587
Electron Energy=-60.0338
Total Energy=-40.5751

Suv_RA
        0  0.346262 -0.346262         0
        0         0         0         0
        0         0         0         0

coef Xuv
-18 -18
-18 -18


gradient (gamma part)
       0  5.41613 -5.41613        0
       0        0        0        0
       0        0        0        0

Yab
-1.5 -1.5
-1.5 -1.5

e term 1
-6.23272  6.23272
       0        0
       0        0

e term 2
-8.12419  8.12419
       0        0
       0        0

e grad total
-14.3569  14.3569
       0        0
       0        0


gradient (Nuclear part)
 13.915 -13.915
      0       0
      0       0

TOTAL GRAD
-0.441958  0.441958
        0         0
        0         0
```

The derivations and figuring out all the indexing took me a long time, and I got some support from the other students, so the code is especially messy, but it does give output that matches the sample output for `H2`. I'll discuss how I got to the outcomes below and attach the pages of hand-written math at the end!

I found that the coefficient to the overlap deriviative `Xuv` is `Ptot_uv * (Ba + Bb)` which ended up not being too hard to figure out since only terms where the AOs are on different atoms are non-zero and these are the only terms that multiply `Suv`. 

The derivative of the overlap wrt coordinates was much harder to understand for me until I realized that the integrals shown in the lab slides are the exact same form as our overlap integral. Once this clicked, it was just a matter of doing all the book-keeping correct to make sure I was indexing the right AOs, primatives, and atoms. My math is pretty weak and I am much more comfortable thinking in code, so once I made the connection that I just needed to call my `__inner_overlap` function twice with different `la` values that make things a lot simpler. 

The gamma terms were much more complex. Most of my time was spent on the gamma derivative. With the help of the lab slides, I came to a solution for the derivative and implemented it (though my implementation is laughibly inefficient... and the code is all tossed in the main function unfortunately). I found the `Yab` matrix coefficient to gamma derivative `Paa*Pbb - ZB*Paa - ZA*Pbb - SUM(p_a^2 + p_b^2) over all pairs uv`. The final summation term I got right away, but the other terms required some guess and check.

My output for `H2` is fully correct, but I couldn't get all of my terms for `HF` or `HO`. My nuclear gradient is always correct, as is my gamma derivative matrix. I suspect my `Yab` and `Xuv` matrices are also correct, but my overlap derivative matrix is not correct.

For example, here is my overlap derivative for `HF`.
My output:

```
Suv_RA
        0  0.357071         0  0.710104  0.710104 -0.357071         0         0         0         0         0         0         0         0         0 -0.710104         0         0         0         0 -0.710104         0         0         0         0
        0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0
        0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0
```

Correct output:
```
Suv_RA (A is the center of u)
        0   0.3571  -0.1807        0        0  -0.3571        0        0        0        0   0.1807        0        0        0        0        0        0        0        0        0        0        0        0        0        0
        0        0        0   0.2044        0        0        0        0        0        0        0        0        0        0        0  -0.2044        0        0        0        0        0        0        0        0        0
        0        0        0        0   0.2044        0        0        0        0        0        0        0        0        0        0        0        0        0        0        0  -0.2044        0        0        0        0
```

The element at `(0,1)` and `(0, 5)` are fine (I'm pretty sure its because they are just H related elements). But clearly I'm not treating the p orbitals fully correctly since I'm missing/incorrect with the other matrix elements.

Here is my full output for `HF`:

```
Nuclear Repulsion Energy=110.742
Electron Energy=-873.141
Total Energy=-762.399

Suv_RA
        0  0.357071         0  0.710104  0.710104 -0.357071         0         0         0         0         0         0         0
0         0 -0.710104         0         0         0         0 -0.710104         0         0         0         0
        0         0         0         0         0         0         0         0         0         0         0         0         0
0         0         0         0         0         0         0         0         0         0         0         0
        0         0         0         0         0         0         0         0         0         0         0         0         0
0         0         0         0         0         0         0         0         0         0         0         0

coef Xuv
-13.9636 -19.4684  42.5338       -0       -0
-19.4684 -145.519 -22.8986       -0       -0
 42.5338 -22.8986 -105.972       -0       -0
      -0       -0       -0     -156       -0
      -0       -0       -0       -0     -156


gradient (gamma part)
       0  5.76911 -5.76911        0
       0        0        0        0
       0        0        0        0

Yab
-1.25061 -7.52514
-7.52514 -55.6991

e term 1
-6.95159  6.95159
       0        0
       0        0

e term 2
-43.4133  43.4133
       0        0
       0        0

e grad total
-50.3649  50.3649
       0        0
       0        0


gradient (Nuclear part)
 64.3851 -64.3851
       0        0
       0        0

TOTAL GRAD
 14.0202 -14.0202
       0        0
       0        0
```

The missing terms in my overlap derivative make the electronic component of the gradient smaller than it should be (`-50` instead of `-58` for the `(0,0)` term).

### Derivations

![overlap](overlap_term.png)

![Yab](Yab.png)

![gamma](gamma_term.png)

