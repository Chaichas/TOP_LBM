/********************  HEADERS  *********************/
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include "lbm_comm.h"

/*******************  FUNCTION  *********************/
int lbm_helper_pgcd(int a, int b)
{
	int c;
	while(b!=0)
	{
		c = a % b;
		a = b;
		b = c;
	}
	return a;
}

/*******************  FUNCTION  *********************/
/**
 * Affiche la configuation du lbm_comm pour un rank donné
 * @param mesh_comm Configuration à afficher
**/
void  lbm_comm_print( lbm_comm_t *mesh_comm )
{
	int rank ;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	printf( " RANK %d ( LEFT %d RIGHT %d TOP %d BOTTOM %d CORNER %d, %d, %d, %d ) ( POSITION %d %d ) (WH %d %d ) \n", rank,
									    mesh_comm->left_id,
									    mesh_comm->right_id,
										mesh_comm->top_id,
									    mesh_comm->bottom_id,
										mesh_comm->corner_id[0],
		 								mesh_comm->corner_id[1],
		 								mesh_comm->corner_id[2],
		 		 						mesh_comm->corner_id[3],
									    mesh_comm->x,
									    mesh_comm->y,
									    mesh_comm->width,
									    mesh_comm->height );
}

/*******************  FUNCTION  *********************/
int helper_get_rank_id(int nb_x,int nb_y,int rank_x,int rank_y)
{
	if (rank_x < 0 || rank_x >= nb_x)
		return -1;
	else if (rank_y < 0 || rank_y >= nb_y)
		return -1;
	else
		return (rank_x + rank_y * nb_x);
}

/*******************  FUNCTION  *********************/
/**
 * Initialise un lbm_comm :
 * - Voisins
 * - Taille du maillage local
 * - Position relative
 * @param mesh_comm MeshComm à initialiser
 * @param rank Rank demandant l'initalisation
 * @param comm_size Taille totale du communicateur
 * @param width largeur du maillage
 * @param height hauteur du maillage
**/
void lbm_comm_init( lbm_comm_t * mesh_comm, int rank, int comm_size, int width, int height )
{
	//vars
	int nb_x;
	int nb_y;
	int rank_x;
	int rank_y;
	//Version incoherente
	/* nb_y = lbm_helper_pgcd(comm_size,width);
	nb_x = comm_size/nb_y; */
	
	//Version coherente (horizontal)
	/* nb_y = lbm_helper_pgcd(comm_size,height);
	nb_x = comm_size/nb_y;*/
	
	// Version verticale (decoupage)
	nb_x = lbm_helper_pgcd(comm_size,width);
	nb_y = comm_size / nb_x;

	//Horizontal
	// nb_y = lbm_helper_pgcd(comm_size,height);
	// nb_x = comm_size / nb_x;

	//check
	assert(nb_x * nb_y == comm_size); //valid for both verticla and horizontal decomposition
	if (height % nb_y != 0) /*this check is only valid for vertical decomposition
		to be modified if we have an horizontal decomposition*/
		fatal("Can't get a 2D cut for current problem size and number of processes.");

	//calc current rank position (ID)
	rank_x = rank % nb_x; //column id
	rank_y = rank / nb_x; //line id

	//setup nb
	mesh_comm->nb_x = nb_x;
	mesh_comm->nb_y = nb_y;

	//setup size (+2 for ghost cells on border)
	mesh_comm->width = width / nb_x + 2;
	mesh_comm->height = height / nb_y + 2;

	//setup position
	mesh_comm->x = rank_x * width / nb_x;
	mesh_comm->y = rank_y * height / nb_y;
	
	// Compute neighbour nodes id
	mesh_comm->left_id  = helper_get_rank_id(nb_x,nb_y,rank_x - 1,rank_y);
	mesh_comm->right_id = helper_get_rank_id(nb_x,nb_y,rank_x + 1,rank_y);
	mesh_comm->top_id = helper_get_rank_id(nb_x,nb_y,rank_x,rank_y - 1);
	mesh_comm->bottom_id = helper_get_rank_id(nb_x,nb_y,rank_x,rank_y + 1);
	mesh_comm->corner_id[CORNER_TOP_LEFT] = helper_get_rank_id(nb_x,nb_y,rank_x - 1,rank_y - 1);
	mesh_comm->corner_id[CORNER_TOP_RIGHT] = helper_get_rank_id(nb_x,nb_y,rank_x + 1,rank_y - 1);
	mesh_comm->corner_id[CORNER_BOTTOM_LEFT] = helper_get_rank_id(nb_x,nb_y,rank_x - 1,rank_y + 1);
	mesh_comm->corner_id[CORNER_BOTTOM_RIGHT] = helper_get_rank_id(nb_x,nb_y,rank_x + 1,rank_y + 1);

	//To be checked later
	//if more than 1 on y, need transmission buffer
	if (nb_y > 1)
	{
		mesh_comm->buffer_upper_send = malloc(sizeof(double) * DIRECTIONS * width / nb_x);
		mesh_comm->buffer_upper_recv = malloc(sizeof(double) * DIRECTIONS * width / nb_x);
		mesh_comm->buffer_lower_send = malloc(sizeof(double) * DIRECTIONS * width / nb_x);
		mesh_comm->buffer_lower_recv = malloc(sizeof(double) * DIRECTIONS * width / nb_x);
	} else {
		mesh_comm->buffer_upper_send = NULL;
		mesh_comm->buffer_upper_recv = NULL;
		mesh_comm->buffer_lower_send = NULL;
		mesh_comm->buffer_lower_recv = NULL;
	}

	//if debug print comm
	#ifndef NDEBUG
	lbm_comm_print( mesh_comm );
	#endif
}


/*******************  FUNCTION  *********************/
/**
 * Libere un lbm_comm
 * @param mesh_comm MeshComm à liberer
**/
void lbm_comm_release( lbm_comm_t * mesh_comm )
{
	mesh_comm->x = 0;
	mesh_comm->y = 0;
	mesh_comm->width = 0;
	mesh_comm->height = 0;
	mesh_comm->right_id = -1;
	mesh_comm->left_id = -1;
	if (mesh_comm->buffer_upper_send != NULL)
		free(mesh_comm->buffer_upper_send);
	mesh_comm->buffer_upper_send = NULL;
	if (mesh_comm->buffer_upper_recv != NULL)
		free(mesh_comm->buffer_upper_recv);
	mesh_comm->buffer_upper_recv = NULL;
	if (mesh_comm->buffer_lower_send != NULL)
		free(mesh_comm->buffer_lower_send);
	mesh_comm->buffer_lower_send = NULL;
	if (mesh_comm->buffer_lower_recv != NULL)
		free(mesh_comm->buffer_lower_recv);
	mesh_comm->buffer_lower_recv = NULL;
}

/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh_to_process Mesh a utiliser lors de l'échange des mailles fantomes
**/
/* void lbm_comm_sync_ghosts_horizontal( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int x )
{
	//vars
	MPI_Status status;

	//if target is -1, no comm
	if (target_rank == -1)
		return;

	//int y;

	switch (comm_type)
	{
		case COMM_SEND:
			//for( y = 0 ; y < mesh->height-2 ; y++ )
				//MPI_Send( &Mesh_get_col( mesh_to_process, x )[y], DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
			MPI_Send( Mesh_get_col( mesh_to_process, x ), DIRECTIONS*(mesh->height-2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
			break;
		case COMM_RECV:
			//for( y = 0 ; y < mesh->height-2 ; y++ )
				//MPI_Recv(  &Mesh_get_col( mesh_to_process, x )[y], DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,&status);
			MPI_Recv( Mesh_get_col( mesh_to_process, x ),  DIRECTIONS*(mesh->height-2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,&status);
			break;
		default:
			fatal("Unknown type of communication.");
	}
}
*/
void lbm_comm_async_ghosts_horizontal( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int x, int Id_req)
{
	//if target is -1, no comm
	if (target_rank == -1)
		return;

	switch (comm_type)
	{
		case COMM_SEND:
			MPI_Isend( Mesh_get_col( mesh_to_process, x ), DIRECTIONS*(mesh->height-2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &mesh->requests[Id_req]);
			break;
		
		case COMM_RECV:
			MPI_Irecv( Mesh_get_col( mesh_to_process, x ),  DIRECTIONS*(mesh->height-2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &mesh->requests[Id_req]);
			break;
		default:
			fatal("Unknown type of communication.");
	}
}


/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh_to_process Mesh a utiliser lors de l'échange des mailles fantomes
**/
/* void lbm_comm_sync_ghosts_diagonal( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int x ,int y)
{
	//vars
	MPI_Status status;

	//if target is -1, no comm
	if (target_rank == -1)
		return;

	switch (comm_type)
	{
		case COMM_SEND:
			MPI_Send( Mesh_get_cell( mesh_to_process, x, y ), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
			break;
		case COMM_RECV:
			MPI_Recv( Mesh_get_cell( mesh_to_process, x, y ), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &status);
			break;
		default:
			fatal("Unknown type of communication.");
	}
} */

void lbm_comm_async_ghosts_diagonal( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int x ,int y, int Id_req)
{
	//if target is -1, no comm
	if (target_rank == -1)
		return;

	switch (comm_type)
	{
		case COMM_SEND:
			MPI_Isend( Mesh_get_cell( mesh_to_process, x, y ), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &mesh->requests[Id_req]);
			break;
		case COMM_RECV:
			MPI_Irecv( Mesh_get_cell( mesh_to_process, x, y ), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &mesh->requests[Id_req]);
			break;
		default:
			fatal("Unknown type of communication.");
	}
}

/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh_to_process Mesh a utiliser lors de l'échange des mailles fantomes
**/
/*void lbm_comm_sync_ghosts_vertical( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int y)
{
	//vars
	MPI_Status status;
	int x;
	//int k;

	//if target is -1, no comm
	if (target_rank == -1)
		return;

	if(comm_type==COMM_SEND) {
			
			for ( x = 1 ; x < mesh_to_process->width - 1 ; x++)
				memcpy(&mesh->buffer[(x-1)*DIRECTIONS], Mesh_get_cell(mesh_to_process, x, y), sizeof(double) * DIRECTIONS);
			MPI_Send( &mesh->buffer[0], DIRECTIONS * (mesh_to_process->width - 2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
				//for ( k = 0 ; k < DIRECTIONS ; k++)
					//MPI_Send( &Mesh_get_cell(mesh_to_process, x, y)[k], 1, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
	}			
	else if(comm_type==COMM_RECV) {

			MPI_Recv( &mesh->buffer[0], DIRECTIONS * (mesh_to_process->width - 2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,&status);
			
			for ( x = 1 ; x < mesh_to_process->width - 1 ; x++)
				memcpy(&Mesh_get_cell(mesh_to_process, x, y)[0], &mesh->buffer[(x-1)*DIRECTIONS], sizeof(double) * DIRECTIONS);
				//for ( k = 0 ; k < DIRECTIONS ; k++)
				//	MPI_Recv( &Mesh_get_cell(mesh_to_process, x, y)[k], DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,&status);
			
	}
	else {
			fatal("Unknown type of communication.");
	}
} 
*/
void lbm_comm_async_ghosts_vertical( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int y, int Id_req)
{
	//vars
	//MPI_Status status;
	int x;

	//if target is -1, no comm
	if (target_rank == -1)
		return;
	if(comm_type==COMM_SEND) {
			if(Id_req%8 == TOP_LINE) {
				for ( x = 1 ; x < mesh_to_process->width - 1 ; x++)
					memcpy(&mesh->buffer_upper_send[(x-1)*DIRECTIONS], Mesh_get_cell(mesh_to_process, x, y), sizeof(double) * DIRECTIONS);
				MPI_Isend( &mesh->buffer_upper_send[0], DIRECTIONS * (mesh_to_process->width - 2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &mesh->requests[Id_req]);
			}
			else if(Id_req%8 == BOTTOM_LINE) {
				for ( x = 1 ; x < mesh_to_process->width - 1 ; x++)
					memcpy(&mesh->buffer_lower_send[(x-1)*DIRECTIONS], Mesh_get_cell(mesh_to_process, x, y), sizeof(double) * DIRECTIONS);
				MPI_Isend( &mesh->buffer_lower_send[0], DIRECTIONS * (mesh_to_process->width - 2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &mesh->requests[Id_req]);
			}
			else {
				fatal("Wrong  bla bla");
			}
	}			
	else if(comm_type==COMM_RECV) {
			if(Id_req%8 == TOP_LINE) {
				MPI_Irecv( &mesh->buffer_upper_recv[0], DIRECTIONS * (mesh_to_process->width - 2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &mesh->requests[Id_req]);
				for ( x = 1 ; x < mesh_to_process->width - 1 ; x++)
					memcpy(&Mesh_get_cell(mesh_to_process, x, y)[0], &mesh->buffer_upper_recv[(x-1)*DIRECTIONS], sizeof(double) * DIRECTIONS);
			} else if(Id_req%8 == BOTTOM_LINE) {
				MPI_Irecv( &mesh->buffer_lower_recv[0], DIRECTIONS * (mesh_to_process->width - 2), MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &mesh->requests[Id_req]);
				for ( x = 1 ; x < mesh_to_process->width - 1 ; x++)
					memcpy(&Mesh_get_cell(mesh_to_process, x, y)[0], &mesh->buffer_lower_recv[(x-1)*DIRECTIONS], sizeof(double) * DIRECTIONS);
			}
			else {
				fatal("Wrong  bla bla");
			}
	}
	else {
			fatal("Unknown type of communication.");
	}
} 

/*******************  FUNCTION  *********************/
/* void lbm_comm_ghost_exchange(lbm_comm_t * mesh, Mesh *mesh_to_process )
{
	//vars
	int rank;

	//get rank
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (rank % 2 == 0){
		//Left to right phase : on reçoit à droite et on envoie depuis la gauche
		lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->right_id,mesh->width - 2);
		lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->left_id,0);

		//prevend comm mixing to avoid bugs
		MPI_Barrier(MPI_COMM_WORLD);
	
		// Right to left phase : on reçoit à gauche et on envoie depuis la droite
		lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->left_id,1);
		lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->right_id,mesh->width - 1);

		//prevend comm mixing to avoid bugs
		MPI_Barrier(MPI_COMM_WORLD);
		
		//top to bottom : on reçoit en bas et on envoie depuis le hauteur
		lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->bottom_id,mesh->height - 2);
		lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->top_id,0);

		//prevend comm mixing to avoid bugs
		MPI_Barrier(MPI_COMM_WORLD);

		// Right to left phase : on reçoit en haut et on envoie depuis le bas
		lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->top_id,1);
		lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->bottom_id,mesh->height - 1);

		//prevend comm mixing to avoid bugs
		MPI_Barrier(MPI_COMM_WORLD);

		//top left
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_LEFT],1,1);
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 1,mesh->height - 1);

		//prevend comm mixing to avoid bugs
		MPI_Barrier(MPI_COMM_WORLD);

		//bottom left
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_LEFT],1,mesh->height - 2);
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 1,0);

		//prevend comm mixing to avoid bugs
		MPI_Barrier(MPI_COMM_WORLD);

		//top right
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 2,1);
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_LEFT],0,mesh->height - 1);

		//prevend comm mixing to avoid bugs
		MPI_Barrier(MPI_COMM_WORLD);

		//bottom left
		//lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_LEFT],1,mesh->height - 2);
		//lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 1,0);

		//prevend comm mixing to avoid bugs
		//MPI_Barrier(MPI_COMM_WORLD);

		//bottom right
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 2,mesh->height - 2);
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_LEFT],0,0);

		//prevend comm mixing to avoid bugs
		MPI_Barrier(MPI_COMM_WORLD);

		// Right to left phase : on reçoit à gauche et on envoie depuis la droite
		//lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->left_id,1);
		//lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->right_id,mesh->width - 1);
		
		//wait for IO to finish, VERY important, do not remove.
		//FLUSH_INOUT();

	} else{
		
		lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->left_id,0);
		lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->right_id,mesh->width - 2);		

		MPI_Barrier(MPI_COMM_WORLD);
		
		lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->right_id,mesh->width - 1);
		lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->left_id,1);

		
		MPI_Barrier(MPI_COMM_WORLD);
		
		lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->top_id,0);
		lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->bottom_id,mesh->height - 2);
	

		MPI_Barrier(MPI_COMM_WORLD);

		lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->bottom_id,mesh->height - 1);
		lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->top_id,1);
		

		MPI_Barrier(MPI_COMM_WORLD);

		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 1,mesh->height - 1);
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_LEFT],1,1);

		
		MPI_Barrier(MPI_COMM_WORLD);

		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 1,0);
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_LEFT],1,mesh->height - 2);		

		
		MPI_Barrier(MPI_COMM_WORLD);

		
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_LEFT],0,mesh->height - 1);
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 2,1);
		

		MPI_Barrier(MPI_COMM_WORLD);

		
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_LEFT],0,0);
		lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 2,mesh->height - 2);

		
		MPI_Barrier(MPI_COMM_WORLD);
		
	}
	
} */

void lbm_comm_async_ghost_exchange(lbm_comm_t * mesh, Mesh *mesh_to_process )
{
	//vars
	int rank;

	//get rank
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (rank % 2 == 0){
		
		//top to bottom : on reçoit en bas et on envoie depuis le hauteur
		lbm_comm_async_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->bottom_id,mesh->height - 2, TOP_LINE);
		lbm_comm_async_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->top_id,0, BOTTOM_LINE*2);

		// Right to left phase : on reçoit en haut et on envoie depuis le bas
		lbm_comm_async_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->top_id,1, BOTTOM_LINE);
		lbm_comm_async_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->bottom_id,mesh->height - 1, TOP_LINE*2);

		//Left to right phase : on reçoit à droite et on envoie depuis la gauche
		lbm_comm_async_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->right_id,mesh->width - 2, LEFT_COL);
		lbm_comm_async_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->left_id,0, RIGHT_COL*2);
	
		// Right to left phase : on reçoit à gauche et on envoie depuis la droite
		lbm_comm_async_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->left_id,1, RIGHT_COL);
		lbm_comm_async_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->right_id,mesh->width - 1, LEFT_COL*2);

		//top left
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_LEFT],1,1, CORNER_BOTTOM_RIGHT);
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 1,mesh->height - 1, 8);

		//bottom left
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_LEFT],1,mesh->height - 2, CORNER_TOP_RIGHT);
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 1,0, CORNER_BOTTOM_LEFT*2);

		//top right
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 2,1, CORNER_BOTTOM_LEFT);
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_LEFT],0,mesh->height - 1, CORNER_TOP_RIGHT*2);

		//bottom right
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 2,mesh->height - 2, CORNER_TOP_LEFT);
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_LEFT],0,0, CORNER_BOTTOM_RIGHT*2);

	} else{
		
		lbm_comm_async_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->left_id,0, RIGHT_COL*2);
		lbm_comm_async_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->right_id,mesh->width - 2, LEFT_COL);		
		
		lbm_comm_async_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->right_id,mesh->width - 1,LEFT_COL*2);
		lbm_comm_async_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->left_id,1,RIGHT_COL);
		
		lbm_comm_async_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->top_id,0,BOTTOM_LINE*2);
		lbm_comm_async_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->bottom_id,mesh->height - 2, TOP_LINE);

		lbm_comm_async_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->bottom_id,mesh->height - 1,TOP_LINE*2);
		lbm_comm_async_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->top_id,1,BOTTOM_LINE);

		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 1,mesh->height - 1,8);
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_LEFT],1,1,CORNER_BOTTOM_RIGHT);

		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 1,0,CORNER_BOTTOM_LEFT*2);
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_LEFT],1,mesh->height - 2,CORNER_TOP_RIGHT);		

		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_LEFT],0,mesh->height - 1,CORNER_TOP_RIGHT*2);
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 2,1,CORNER_BOTTOM_LEFT);
		
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_LEFT],0,0,CORNER_BOTTOM_RIGHT*2);
		lbm_comm_async_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 2,mesh->height - 2,CORNER_TOP_LEFT);
		
	}
	MPI_Waitall(16, mesh->requests , MPI_STATUSES_IGNORE);
	
}

/*******************  FUNCTION  *********************/
/**
 * Rendu du mesh en effectuant une réduction a 0
 * @param mesh_comm MeshComm à utiliser
 * @param temp Mesh a utiliser pour stocker les segments
**/
void save_frame_all_domain( FILE * fp, Mesh *source_mesh, Mesh *temp )
{
	//vars
	int i = 0;
	int comm_size, rank ;
	MPI_Status status;

	//get rank and comm size
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	/* If whe have more than one process */
	if( 1 < comm_size )
	{
		if( rank == 0 )
		{
			/* Rank 0 renders its local Mesh */
			save_frame(fp,source_mesh);
			/* Rank 0 receives & render other processes meshes */
			for( i = 1 ; i < comm_size ; i++ )
			{
				MPI_Recv( temp->cells, source_mesh->width  * source_mesh->height * DIRECTIONS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
				save_frame(fp,temp);
			}
		} else {
			/* All other ranks send their local mesh */
			MPI_Send( source_mesh->cells, source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		}
	} else {
		/* Only 0 renders its local mesh */
		save_frame(fp,source_mesh);
	}

}

