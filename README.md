# filterSqp
Implementation of the FilterSQP constrained optimization algorithm. Using the FADBAD++ library for automatic 
differentiation.

2015/Sep - This is still very much WIP, and there is nothing related to constrained optimization yet. Just implementing trust-region optimization first, adapting my old Python code from the [Cosisco project](https://github.com/nlw0/corisco).

But we do have preliminary results. Here are some basic demos of the solution of classic functions.

![Himmelblau function optimization from different starting points, using a fixed step length trust region method](https://raw.githubusercontent.com/nlw0/filterSqp/master/demo/himmelblau_min.png)

On the Rosenbrock "banana" function":
![Gradient descent, smooth but slow, too many iterations on the valley.](https://github.com/nlw0/filterSqp/blob/master/demo/gradient_descent.png)

![Newton's method. Faster but clumsy, better to have a smoother track.](https://github.com/nlw0/filterSqp/blob/master/demo/newtonA.png)

![Fixed step trust region with same starting point from Newton's method. Much smoother, following the valley.](https://github.com/nlw0/filterSqp/blob/master/demo/trust_simple_A.png)

![Fixed step trust region with same starting point from gradient descent. Very similar track, but far fewer iterations.](https://github.com/nlw0/filterSqp/blob/master/demo/trust_simple_B.png)
