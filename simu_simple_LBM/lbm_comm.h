#ifndef MESH_COMM_H
#define MESH_COMM_H

/********************  HEADERS  *********************/
#include <mpi.h>
#include <stdlib.h>
#include "lbm_struct.h"

/*******************  DEFINITIONS  ******************/
/** Definition de l'ID du processus maître. **/
#define RANK_MASTER 0

/** Definition pour communication asynchrone verticale. **/

/*********************  ENUM  ***********************/
/**
 * Definition des différents type de cellule pour savoir quel traitement y appliquer
 * lors du calcul.
**/
/* typedef enum lbm_corner_pos_e
{
	CORNER_TOP_LEFT = 0,
	CORNER_TOP_RIGHT = 1,
	CORNER_BOTTOM_LEFT = 2,
	CORNER_BOTTOM_RIGHT = 3,
} lbm_corner_pos_t; */

typedef enum lbm_pos_e
{
	CORNER_TOP_LEFT = 0,
	CORNER_TOP_RIGHT = 1,
	CORNER_BOTTOM_LEFT = 2,
	CORNER_BOTTOM_RIGHT = 3,
	LEFT_COL = 4,
	RIGHT_COL = 5,
	TOP_LINE = 6,
	BOTTOM_LINE = 7,
} lbm_pos_t;

/*********************  ENUM  ***********************/
typedef enum lbm_comm_type_e
{
	COMM_SEND,
	COMM_RECV
} lbm_comm_type_t;

/********************  STRUCT  **********************/
/**
 * Structure utilisée pour stoquer les informations relatives aux communications.
**/
typedef struct lbm_comm_t_s
{
	/** Position de la maille locale dans le maillage global (origine). **/
	int x;
	int y;
	/** Taille de la maille locale. **/
	int width;
	int height;
	int nb_x;
	int nb_y;
	/** Id du voisin de droite, -1 si aucun. **/
	int right_id;
	/** Id du voisin de gauche, -1 si aucun. **/
	int left_id;
	int top_id;
	int bottom_id;
	int corner_id[4];
	/** Requète asynchrone en cours. **/
	//MPI_Request requests[32];
	MPI_Request requests[16];
	int id_req;
	//lbm_mesh_cell_t buffer;
	lbm_mesh_cell_t buffer_upper_send;
	lbm_mesh_cell_t buffer_upper_recv;
	lbm_mesh_cell_t buffer_lower_send;
	lbm_mesh_cell_t buffer_lower_recv;
} lbm_comm_t;

/*******************  FUNCTION  *********************/
static inline int lbm_comm_width( lbm_comm_t *mc )
{
	return mc->width;
}

/*******************  FUNCTION  *********************/
static inline int lbm_comm_height( lbm_comm_t *mc )
{
	return mc->height;
}

/*******************  FUNCTION  *********************/
void lbm_comm_init( lbm_comm_t * mesh, int rank, int comm_size, int width, int height );
void lbm_comm_release( lbm_comm_t * mesh );
void  lbm_comm_print( lbm_comm_t *mesh );

/*******************  FUNCTION  *********************/
void lbm_comm_sync_ghosts_wait( lbm_comm_t * mesh );
//void lbm_comm_ghost_exchange(lbm_comm_t * mesh, Mesh *mesh_to_process );
void lbm_comm_async_ghost_exchange(lbm_comm_t * mesh, Mesh *mesh_to_process );

/*******************  FUNCTION  *********************/
void save_frame_all_domain( FILE * fp, Mesh *source_mesh, Mesh *temp );

#endif
