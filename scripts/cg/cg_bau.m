% x = cg_bau(A, b, x, maxit=0, f=[])
function x = cg_bau(A, b, x, maxit=0, f=[])


[m,n] = size(A);
r = b - A * x;
p = r;

if maxit==0
	maxit=n;
end

i=0;
tol = sqrt( transpose(r) * r );

while i<maxit && tol>1e-6 
	t = A * p;

	u = transpose(r) * r;
	alpha = u / ( transpose(p) * t );
	x += alpha * p;
	r -= alpha * t;
	tol = transpose(r) * r;
	beta = tol / u;
	p = r + beta *p;

	tol = sqrt(tol);
	printf("%i %16.20e %16.20e %16.20e\n", i, tol, alpha, beta);
	
	%anorm = sqrt(transpose(x - f) * A * (x-f));
	i = i+1;
end
