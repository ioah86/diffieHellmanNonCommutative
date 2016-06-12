Diffie Hellman -- Non commutative version
==================================================

In this repository, one can find a version of the Diffie-Hellman
protocol using a non-commutative domain. A publication about this
protocol is being written up and will appear soon.

In order to compile and run the protocol as stated in the paper
"A Diffie-Hellman-like Key Exchange Protocol Based on Multivariate Ore
Polynomials" by Reinhold Burger and Albert Heinle
(http://arxiv.org/abs/1407.1270),
simply use the command

```
make
```

After successful compilation, there will be two new executables in
the folder: `impl` and `tests`. You can run `tests` in order to see if
your compiled version of the code executes the tests correctly. If
this is the case, you know that your version (or your alteration of
the code) passes selected tests. By running `impl`, a simulation of
the protocol between two parties *A* and *B* is run. 

You can use the `ore_algebra` library we created for your own
projects. Simply call

```
make ore_algebra.o
```
to create a `C`-library, and link it to your project. The header file
`ore_algebra.h` is included in the source files.

We are looking forward to feedback from users.
