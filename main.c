#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void apply_filters(char *images_file_name, int *neighbors, int rank, int nr_ps, char *statistics_file_name);

void read_topology(char *filename, int *neighbors, int rank, int nr_ps)
{
	FILE *topologyFile;
	char line[100];
	char separater;
	int line_nr, i, neighbor_rank;
	topologyFile = fopen(filename, "r");
	if (topologyFile == NULL)
	{
		printf("Error opening file!\n");
		exit(0);
	}
	for (i = 0; i < rank; i++)
	{
		if(fgets(line, (nr_ps * 2 + 1), topologyFile) == NULL)
		{
			printf("Error occured while reading file!\n");
			exit(0);
		}
	}
	fscanf(topologyFile, "%d", &line_nr);
	if (rank != line_nr) {
		printf("Wrong line! Check file format.\n");
		exit(0);
	}
	fscanf(topologyFile, "%c", &separater);
	while (separater != '\n')
	{
		fscanf(topologyFile, "%d", &neighbor_rank);
		neighbors[neighbor_rank] = 1;
		if (feof(topologyFile)) {
			break;
		}
		fscanf(topologyFile, "%c", &separater);
	}

	fclose(topologyFile);
}



int main(int argc, char *argv[])
{
	int rank, nr_ps;
	int *neighbors;
	int i;

	MPI_Status status;
	FILE *imagesFile, *staticticsFile;
	

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nr_ps);

	neighbors =(int*) malloc(nr_ps * sizeof(int));
	for (i = 0; i < nr_ps; i++)
	{
		neighbors[i] = 0;
	}
	read_topology(argv[1], neighbors, rank, nr_ps);
	apply_filters(argv[2], neighbors, rank, nr_ps, argv[3]);

	MPI_Finalize();
	return 0;
}