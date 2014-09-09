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
gcc impl.c
```

If you want to do some operations with the underlying ore_algebra, use
the functions as provided in
```
ore_algebra.c
```
for your own project.

We are looking forward to feedback from users.
