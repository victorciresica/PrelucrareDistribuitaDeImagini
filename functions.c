#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SOBEL 7
#define MEAN_REMOVAL 8
#define DONE 13

void image_header(FILE *current_image, FILE *new_image, int *width, int *height)
{
	
	char line[100];
	int i, max_val;

	
	for (i = 0; i < 2; i++)
	{
		if(fgets(line, 100, current_image) == NULL)
		{
			printf("Error writing to file!\n");
			exit(0);
		}
		fprintf(new_image, "%s", line);
	}
	fscanf(current_image, "%d %d\n%d", width, height, &max_val);
	fprintf(new_image, "%d %d\n%d\n", *width, *height, max_val);

}

void scatter_pixels(int *neighbors, int nr_ps, int height, int width,
	int **buffer, int tag, int parent, int *top_edge, int *bottom_edge, int nr_children)
{
	int i, j, count = 0, buf_dim, bonus;
	int start;
	
	buf_dim = height / nr_children;
	bonus = height % nr_children;

	if (buf_dim > 0)
	{
		for (i = 0; i < nr_ps; i++)
		{
			if (neighbors[i] == 1 && i != parent)
			{
				
				start = count * buf_dim;
				if (count == nr_children - 1)
				{
					buf_dim += bonus;
				}
				MPI_Send(&buf_dim, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
				MPI_Send(&width, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
				for (j = 0; j < buf_dim; j++)
				{
					MPI_Send(buffer[start + j], width, MPI_INT, i, tag, MPI_COMM_WORLD);
				}
				// TOP EDGE
				if (count == 0)
					MPI_Send(top_edge, width, MPI_INT, i, tag, MPI_COMM_WORLD);
				else
					MPI_Send(buffer[start - 1], width, MPI_INT, i, tag, MPI_COMM_WORLD);

				// BOTTOM EDGE
				if (count == nr_children - 1)
					MPI_Send(bottom_edge, width, MPI_INT, i, tag, MPI_COMM_WORLD);
				else
				{
					MPI_Send(buffer[start + buf_dim], width, MPI_INT, i, tag, MPI_COMM_WORLD);
				}
				count++;
			}			
		}
	}
}

void filter_pixels(int height, int width, int **buffer, int tag, int *top_edge, int *bottom_edge)
{
	int i, j, sum;
	int **new_buffer;

	new_buffer = (int **) malloc(height * sizeof(int*));
	for (i = 0; i < height; i++)
		new_buffer[i] = (int*) malloc(width * sizeof(int));

	if (tag == SOBEL)
	{
		for (i = 0; i < height; i++)
		{
			for (j = 0; j < width; j++)
			{
				sum = 0;
				if (j != 0)
				{
					sum += 2 * buffer[i][j - 1];
					if (i != 0)
						sum += 1 * buffer[i - 1][j - 1];
					else
						sum += 1 * top_edge[j - 1];
					if (i != height - 1)
						sum += 1 * buffer[i + 1][j - 1];
					else
						sum += 1 * bottom_edge[j - 1];
				}
				if (j != width - 1)
				{
					sum -= 2 * buffer[i][j + 1];
					if (i != 0)
						sum -= buffer[i - 1][j + 1];
					else
						sum -= top_edge[j + 1];
					if (i != height - 1)
						sum -= buffer[i + 1][j + 1];
					else
						sum -= bottom_edge[j + 1];
				}
				sum += 127;
				if (sum < 0)
					sum = 0;
				if (sum > 255)
					sum = 255;
				new_buffer[i][j] = sum;
			}
		}
	}
	if (tag == MEAN_REMOVAL)
	{
		for (i = 0; i < height; i++)
		{
			for (j = 0; j < width; j++)
			{
				sum = 0;
				if (j != 0)
				{
					sum -= buffer[i][j - 1];
					if (i != 0) 
						sum -= buffer[i - 1][j - 1];
					else
						sum -= top_edge[j - 1];
					if (i != height - 1)
						sum -= buffer[i + 1][j - 1];
					else
						sum -= bottom_edge[j - 1];
				}

				if (i != 0)
					sum -= buffer[i - 1][j];
				else
					sum -= top_edge[j];
				if (i != height - 1)
					sum -= buffer[i + 1][j];
				else
					sum -= bottom_edge[j];

				if (j != width - 1)
				{
					sum -= buffer[i][j + 1];
					if (i != 0)
						sum -= buffer[i - 1][j + 1];
					else
						sum -= top_edge[j + 1];
					if (i != height - 1)
						sum -= buffer[i + 1][j + 1];
					else
						sum -= bottom_edge[j + 1];
				}
				sum += 9 * buffer[i][j];
				if (sum < 0)
					sum = 0;
				if (sum > 255)
					sum = 255;
				new_buffer[i][j] = sum;
			}	
		}
	}
	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
			buffer[i][j] = new_buffer[i][j];
		free(new_buffer[i]);
	}
	free(new_buffer);
}

void send_pixels(int** buffer, int height, int width, int parent)
{
	int i;
	for (i = 0; i < height; i++)
	{
		MPI_Send(buffer[i], width, MPI_INT, parent, i, MPI_COMM_WORLD);
	}
}

void gather_pixels(int rank, int nr_ps, int *neighbors, int **buffer,
	int height, int width, int parent, FILE *new_image)
{
	int *order, *lines, *tem_buffer;
	int buf_dim, count = 0;
	int i, j, k;
	MPI_Status status;


	tem_buffer = (int*) malloc(width * sizeof(int));

	order = (int*) malloc(nr_ps * sizeof(int));
	for (i = 0; i < nr_ps; i++)
	{
		if (neighbors[i] == 1 && i != parent)
		{
			order[i] = count;
			count++;
		}
	}

	buf_dim = height / count;

	for (i = 0; i < height; i++)
	{
		MPI_Recv(tem_buffer, width, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		for (j = 0; j < width; j++)
		{
			buffer[order[status.MPI_SOURCE] * buf_dim + status.MPI_TAG][j] = tem_buffer[j];
		}
	}
	if(rank == 0)
	{
		for (i = 0; i < height; i++)
		{
			for (j = 0; j < width; j++)
			{
				fprintf(new_image, "%d\n", buffer[i][j]);
			}
		}
		fclose(new_image);
	}
	else
	{
		send_pixels(buffer, height, width, parent);
	}
}

int scatter_filter_pixels(int rank, int nr_ps, int *neighbors, int width,
	int height, FILE *current_image, FILE *new_image, char *filter, int *tag)
{
	int i, j, buf_dim, lines;
	int **buffer;
	int parent = -1;
	int nr_children = 0;
	int *bottom_edge, *top_edge;
	MPI_Status status;

	if (rank == 0)
	{
		// SET TAG ACCORDING TO THE WANTED FILTER
		if (strcmp(filter, "sobel") == 0)
			*tag = SOBEL;
		else
			if (strcmp(filter, "mean_removal") == 0)
				*tag = MEAN_REMOVAL;
			else
			{
				printf("I can not apply this filter, it's too damn complicated!\n");
				exit(0);
			}

		// ALLOCATE BUFFER AND EDGES
		buffer = (int**) malloc(height * sizeof(int*));
		for (i = 0; i < height; i++)
		{
			buffer[i] = (int*) malloc(width * sizeof(int));
			for (j = 0; j < width; j++)
			{
				fscanf(current_image, "%d", &buffer[i][j]);
			}
		}
		fclose(current_image);
		bottom_edge = (int*) malloc(width * sizeof(int));
		top_edge = (int*) malloc(width * sizeof(int));
		for (j = 0; j < width; j++)
		{
			top_edge[j] = 0;
			bottom_edge[j] = 0;
		}
	}

	// RECEIVE PIXELS
	else
	{
		// RECEIVE THE NUMBER OF LINES
		MPI_Recv(&height, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		parent = status.MPI_SOURCE;
		*tag = status.MPI_TAG;
		if (*tag == DONE)
			return 0;

		MPI_Recv(&width, 1, MPI_INT, parent, *tag, MPI_COMM_WORLD, &status);
		// ALLOCATE BUFFER AND EDGES
		bottom_edge = (int*) malloc(width * sizeof(int));
		top_edge = (int*) malloc(width * sizeof(int));
		buffer = (int**) malloc(height * sizeof(int*));

		for (i = 0; i < height; i++)
		{
			buffer[i] = (int*) malloc(width * sizeof(int));
			MPI_Recv(buffer[i], width, MPI_INT, parent, *tag, MPI_COMM_WORLD, &status);
		}

		MPI_Recv(top_edge, width, MPI_INT, parent, *tag, MPI_COMM_WORLD, &status);
		MPI_Recv(bottom_edge, width, MPI_INT, parent, *tag, MPI_COMM_WORLD, &status);
	}

	// COUNT CHILDREN
	for (i = 0; i < nr_ps; i++)
	{
		if (neighbors[i] == 1 && i != parent)
			nr_children++;
	}
	if (nr_children != 0)
	{
		lines = 0;
		scatter_pixels(neighbors, nr_ps, height, width, buffer, *tag, parent, top_edge, bottom_edge, nr_children);
		gather_pixels(rank, nr_ps, neighbors, buffer, height, width, parent, new_image);
	}
	else
	{
		lines = height;
		filter_pixels(height, width, buffer, *tag, top_edge, bottom_edge);
		send_pixels(buffer, height, width, parent);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for (i = 0; i < height; i++) {
		free(buffer[i]);
	}
	free(buffer);
	if (rank != 0)
	{
		free(top_edge);
		free(bottom_edge);
	}
	return lines;
}

void statistics(int rank, int nr_ps, int *neighbors, int lines, char *statistics_file_name)
{

	int buffer, i, j, parent = -1, nr_children = 0;
	int *statistic_table, *table_buffer;
	FILE *statistics_file;
	MPI_Status status;
	
	statistic_table = (int*) malloc(nr_ps * sizeof(int));
	table_buffer = (int*) malloc(nr_ps * sizeof(int));
	for (i = 0; i < nr_ps; i++) {
		statistic_table[i] = 0;
	}
	statistic_table[rank] = lines;
	if (rank != 0)
	{
		MPI_Recv(&buffer, 1, MPI_INT, MPI_ANY_SOURCE, DONE, MPI_COMM_WORLD, &status);
		parent = status.MPI_SOURCE;
	}
	for (i = 0; i < nr_ps; i++)
		if (neighbors[i] == 1 && i != parent) 
		{
			nr_children++;
			MPI_Send(&rank, 1, MPI_INT, i, DONE, MPI_COMM_WORLD);
			MPI_Send(&rank, 1, MPI_INT, i, DONE, MPI_COMM_WORLD);
		}
	if (nr_children == 0)
	{
		MPI_Send(statistic_table, nr_ps, MPI_INT, parent, DONE, MPI_COMM_WORLD);
	}
	else
	{
		
		for (i = 0; i < nr_children; i++)
		{
			MPI_Recv(table_buffer, nr_ps, MPI_INT, MPI_ANY_SOURCE, DONE, MPI_COMM_WORLD, &status);
			for (j = 0; j < nr_ps; j++)
			{
				if (table_buffer[j] != 0)
					statistic_table[j] = table_buffer[j];
			}
		}
		if (rank != 0)
		{
			MPI_Send(statistic_table, nr_ps, MPI_INT, parent, DONE, MPI_COMM_WORLD);
		}
		else
		{
			statistics_file = fopen(statistics_file_name, "w");
			for (i = 0; i < nr_ps; i++)
				fprintf(statistics_file, "%d: %d\n", i, statistic_table[i]);
			fclose(statistics_file);
		}
	}
	free(statistic_table);
	free(table_buffer);
}

void apply_filters(char *images_file_name, int *neighbors, int rank, int nr_ps, char *statistics_file_name)
{
	char image_name[100], filter[20], new_image_name[100];
	int nr_images, i, image_width = 0, image_height = 0, tag = 3, lines = 0;
	FILE *images_file, *current_image, *new_image;

	if (rank == 0)
	{
		images_file = fopen(images_file_name, "r");
		fscanf(images_file, "%d", &nr_images);

		for (i = 0; i < nr_images; i++)
		{
			fscanf(images_file, "%s %s %s", filter, image_name, new_image_name); 
			current_image = fopen(image_name, "r");
			if(current_image == NULL) 
			{
				printf("Can't open file %s!\n", image_name);
				exit(0);
			}
			new_image = fopen(new_image_name, "w");
			if (new_image == NULL)
			{
				printf("Can't open file %s!\n", image_name);
				exit(0);		
			}
			image_header(current_image, new_image, &image_width, &image_height);
			lines += scatter_filter_pixels(rank, nr_ps, neighbors, image_width, image_height,
			 current_image, new_image, filter, &tag);
		}
		fclose(images_file);
		statistics(rank, nr_ps, neighbors, lines, statistics_file_name);
	}
	else
	{
		while (tag != DONE)
			lines += scatter_filter_pixels(rank, nr_ps, neighbors, image_width, image_height, 
				current_image, new_image, filter, &tag);
		statistics(rank, nr_ps, neighbors, lines, statistics_file_name);
	}
}