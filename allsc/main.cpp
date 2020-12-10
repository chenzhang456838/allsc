#include "Utility.h"
#include "Graph.h"

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("Usage: exe [2]graph-dir\n");
		return 0;
	}

#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#else
	int start, end1, end;
	start = clock();
#endif

	Graph *graph = new Graph(argv[1]);
	graph->read_graph();
	printf("Finished reading graph\n");

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;

#else
	end1 = clock();

#endif

	graph->compute_allSC();

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("IO time: %lld, all_SC ltime: %lld\n", mtime1, mtime);

#else
	end = clock();

	printf("IO time: %lld, all_SC time: %d\n", end1-start, end-start);
#endif

	return 0;
}
