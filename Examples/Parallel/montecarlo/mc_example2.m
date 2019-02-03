# this runs the test on several different configurations and reports timings.
# to see good speedups in parallel, you need to set "reps" to a fairly
# big number (say, 100000 or so).
outfile = "example2.out";
T = 1000;
dim = 5;
reps = 10000;

maxcomputenodes = 2;
timings = zeros(maxcomputenodes,1);
for i = 0:maxcomputenodes;
	tic;
	n_received = montecarlo("tracetest", {T,dim}, reps, outfile, i, 1000);
	t = toc;
	if i > 0
		bar(0:i, n_received/reps);
		xlabel("nodes");
		grid on;
		title("percentage of results received by node: note - the frontend (node 0) does not do calculations");
		drawnow();
		if (i < maxcomputenodes) figure; endif
	endif

	timings(i+1,:) = t;
endfor

printf("\n\nTiming results, second example\n");
for i = 0:maxcomputenodes;
	printf("Number of compute nodes: %d     Time (sec.): %f\n", i, timings(i+1,:));
endfor
