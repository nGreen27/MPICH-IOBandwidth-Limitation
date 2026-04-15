/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#include "mpiimpl.h"
//#include "mpiu_thread.h"
#include "async_io_thread.h"

/*
 *   Global variables
 */

/* Request queue for async_io_thread */
MPIIO_Async_io_queue_t MPIIO_Async_io_queue;


#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE
static MPID_Thread_id_t async_io_thread_id;
static MPID_Thread_mutex_t async_io_mutex;
static MPID_Thread_cond_t async_io_cond_dequeue;
static MPID_Thread_cond_t async_io_cond_enqueue;
#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */

/*
 *   Thread related functions
 */

/* Init queue for async_io_thread */
int MPIIO_Init_async_io_queue(MPIIO_Async_io_queue_t *queue, int queue_size, int elem_size)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Assert(queue != NULL);
#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE

    queue->queue_size = queue_size;
    queue->elem_size = elem_size;
    queue->front = 0;
    queue->rear = 0;
    queue->array = (MPIIO_gen_queue_elem_t *)MPL_malloc (sizeof(MPIIO_gen_queue_elem_t)*queue_size,MPL_MEM_THREAD);
    if (queue->array == NULL) goto fn_fail1;
    for (int i=0; i<queue_size; i++) queue->array[i].data=NULL;
    for (int i=0; i<queue_size; i++) {
        queue->array[i].oper = 0;
        queue->array[i].size = 0;
        queue->array[i].data=(char *)MPL_malloc (elem_size,MPL_MEM_THREAD);
        if (queue->array[i].data == NULL) goto fn_fail2;
    }
fn_exit:
    return mpi_errno;
fn_fail1:
    return MPI_ERR_OTHER;
fn_fail2:
    for (int i=0; i<queue_size; i++) {
        if (queue->array[i].data != NULL) MPL_free(queue->array[i].data);
    }
    MPL_free(queue->array);
    return MPI_ERR_OTHER;

#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */
    return mpi_errno;
}

/* End queue for async_io_thread */
int MPIIO_End_async_io_queue(MPIIO_Async_io_queue_t *queue)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Assert(queue != NULL);
#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE

    if (queue->array != NULL) {
        for (int i=0; i<queue->queue_size; i++) {
            if (queue->array[i].data != NULL) {
                MPL_free(queue->array[i].data);
            }
        }
        MPL_free(queue->array);
    }
    queue->queue_size = 0;
    queue->elem_size = 0;
    queue->front = 0;
    queue->rear = 0;
    
#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */
    return mpi_errno;
}

/* check if queue for async_io_thread is empty */
int MPIIO_elem_size_async_io_queue(MPIIO_Async_io_queue_t *queue)
{
    MPIR_Assert(queue != NULL);
    return (queue->elem_size);
}

/* check if queue for async_io_thread is empty */
int MPIIO_IsEmpty_async_io_queue(MPIIO_Async_io_queue_t *queue)
{
    MPIR_Assert(queue != NULL);
    return (queue->front == queue->rear);
}

/* check if queue for async_io_thread is full */
int MPIIO_IsFull_async_io_queue(MPIIO_Async_io_queue_t *queue)
{
    MPIR_Assert(queue != NULL);
    return ((queue->rear+1)%queue->queue_size == queue->front);
}

/* enqueue element on queue for async_io_thread */
int MPIIO_Enqueue_async_io_queue(MPIIO_Async_io_queue_t *queue, int oper, char *data, int size)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Assert(queue != NULL);
    MPIR_Assert(data != NULL);
    MPIR_Assert(size <= queue->elem_size);

#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE

    /* XXX DJG why is this unlock/lock necessary?  Should we just YIELD here or later?  */
    //MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

    MPID_Thread_mutex_lock(&async_io_mutex, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    while (MPIIO_IsFull_async_io_queue(queue)) {
        PRINTF("This Queue is full\n");
        MPID_Thread_cond_wait(&async_io_cond_dequeue, &async_io_mutex, &mpi_errno);
        MPIR_Assert(!mpi_errno);
    }
    queue->rear = (queue->rear +1)%queue->queue_size;
    PRINTF("enqueued element: pos=%d, op=%d size=%d\n", queue->rear, oper, size);
    /* copy the values directly to the queue element */
    queue->array[queue->rear].oper = oper;
    queue->array[queue->rear].size = size;
    MPIR_Memcpy(queue->array[queue->rear].data, data, size);

    MPID_Thread_mutex_unlock(&async_io_mutex, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    MPID_Thread_cond_signal(&async_io_cond_enqueue, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    //MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */
    return mpi_errno;
}

/* dequeue element from queue for async_io_thread */
/* note: input size -> booked size; ouput size -> data size */
int MPIIO_Dequeue_async_io_queue(MPIIO_Async_io_queue_t *queue, int *oper, char *data, int *size)
{
    int mpi_errno = MPI_SUCCESS;

    MPIR_Assert(queue != NULL);
    MPIR_Assert(oper != NULL);
    MPIR_Assert(data != NULL);
    MPIR_Assert(size != NULL);
    MPIR_Assert((*size) >= queue->elem_size);

#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE

    /* XXX DJG why is this unlock/lock necessary?  Should we just YIELD here or later?  */
    //MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

    MPID_Thread_mutex_lock(&async_io_mutex, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    while (MPIIO_IsEmpty_async_io_queue(queue)) {
        PRINTF("This Queue is empty\n");
        MPID_Thread_cond_wait(&async_io_cond_enqueue, &async_io_mutex, &mpi_errno);
        MPIR_Assert(!mpi_errno);
    }
    queue->front = (queue->front +1)%queue->queue_size;
    /* copy the values directly to the queue element */
    (*oper) = queue->array[queue->front].oper;
    (*size) == queue->array[queue->front].size;
    MPIR_Memcpy(data, queue->array[queue->front].data, (*size));
    PRINTF("dequeued element: pos=%d, op=%d, size=%d\n", queue->front, (*oper), (*size));

    MPID_Thread_mutex_unlock(&async_io_mutex, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    MPID_Thread_cond_signal(&async_io_cond_dequeue, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    //MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */
    
    return mpi_errno;
}

/* wait until queue for async_io_thread is empty */
int MPIIO_WaitEmpty_async_io_queue(MPIIO_Async_io_queue_t *queue)
{
    int mpi_errno = MPI_SUCCESS;

#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE

    /* XXX DJG why is this unlock/lock necessary?  Should we just YIELD here or later?  */
    //MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

    MPID_Thread_mutex_lock(&async_io_mutex, &mpi_errno);

    MPIR_Assert(!mpi_errno);
    while (!MPIIO_IsEmpty_async_io_queue(queue)) {
        PRINTF("This Queue is not empty\n");
        MPID_Thread_cond_wait(&async_io_cond_dequeue, &async_io_mutex, &mpi_errno);
        MPIR_Assert(!mpi_errno);
    }
    PRINTF("This Queue is empty\n");
    MPID_Thread_mutex_unlock(&async_io_mutex, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    //MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */
    
    return mpi_errno;
}

/* main function of the Asynchronous IO thread */
#undef FUNCNAME
#define FUNCNAME async_io_fn
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
static void async_io_fn(void)
{
    int mpi_errno = MPI_SUCCESS;

#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE
    /* Explicitly add CS_ENTER/EXIT since this thread is created from
     * within an internal function and will call NMPI functions
     * directly. */
    //MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    PRINTF("async_io_fn: begin\n");

    /* allocate mem. for the data element */
    int elem_size = MPIIO_elem_size_async_io_queue(&MPIIO_Async_io_queue);
    PRINTF("async_io_fn:MAX elem size = %d\n",elem_size);
    char *data = (char *)MPL_malloc (elem_size,MPL_MEM_THREAD);
    int oper = 0;
    do {
        int size = elem_size;
        mpi_errno = MPIIO_Dequeue_async_io_queue(&MPIIO_Async_io_queue, &oper, (char *)data, &size);
        MPIR_Assert(!mpi_errno);
        PRINTF("async_io_fn: Dequeue oper = %d, size = %d\n",oper, size);
        /* if valid operation, execute it */
        if (oper == 0) {
            break;
        }
        Async_exec_op(oper, (void *)data, size);
    } while (oper != 0);

    /* free mem. for the data element */
    PRINTF("async_io_fn: free temp data\n");
    MPL_free(data);

    PRINTF("async_io_fn: end\n");
    //MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */
    return;
}

/* Init the Asynchronous IO thread */
#undef FUNCNAME
#define FUNCNAME MPIIO_Init_async_io_thread
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIIO_Init_async_io_thread(int queue_size, int elem_size)
{
#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE
    int mpi_errno = MPI_SUCCESS;
    int err = 0;

    PRINTF("MPIIO_Init_async_io_thread: begin\n");

    err = MPIIO_Init_async_io_queue(&MPIIO_Async_io_queue, queue_size, elem_size);
    MPIR_ERR_CHKANDJUMP1(err, mpi_errno, MPI_ERR_OTHER, "**queue_create", "**queue_create %s", strerror(err));

    MPID_Thread_cond_create(&async_io_cond_dequeue, &err);
    MPIR_ERR_CHKANDJUMP1(err, mpi_errno, MPI_ERR_OTHER, "**cond_create", "**cond_create %s", strerror(err));

    MPID_Thread_cond_create(&async_io_cond_enqueue, &err);
    MPIR_ERR_CHKANDJUMP1(err, mpi_errno, MPI_ERR_OTHER, "**cond_create", "**cond_create %s", strerror(err));

    MPID_Thread_mutex_create(&async_io_mutex, &err);
    MPIR_ERR_CHKANDJUMP1(err, mpi_errno, MPI_ERR_OTHER, "**mutex_create", "**mutex_create %s", strerror(err));
    
    MPID_Thread_create((MPID_Thread_func_t) async_io_fn, NULL, &async_io_thread_id, &err);
    MPIR_ERR_CHKANDJUMP1(err, mpi_errno, MPI_ERR_OTHER, "**mutex_create", "**mutex_create %s", strerror(err));
    
    PRINTF("MPIIO_Init_async_io_thread: end\n");

 fn_exit:
    return mpi_errno;
 fn_fail:
    goto fn_exit;
#else
    return MPI_SUCCESS;
#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */
}

/* End the Asynchronous IO thread */
#undef FUNCNAME
#define FUNCNAME MPIO_Finalize_async_io_thread
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIO_Finalize_async_io_thread(void)
{
    int mpi_errno = MPI_SUCCESS;
#if MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE

    PRINTF("MPIO_Finalize_async_io_thread: begin\n");

    int data = 0; /* End operation */
    mpi_errno = MPIIO_Enqueue_async_io_queue(&MPIIO_Async_io_queue, 0, (char *)&data, sizeof(int));
    MPIR_Assert(!mpi_errno);

    mpi_errno = MPIIO_WaitEmpty_async_io_queue(&MPIIO_Async_io_queue);
    MPIR_Assert(!mpi_errno);

    MPID_Thread_cond_destroy(&async_io_cond_dequeue, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    MPID_Thread_cond_destroy(&async_io_cond_enqueue, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    MPID_Thread_mutex_destroy(&async_io_mutex, &mpi_errno);
    MPIR_Assert(!mpi_errno);

    mpi_errno = MPIIO_End_async_io_queue(&MPIIO_Async_io_queue);
    MPIR_Assert(!mpi_errno);

    PRINTF("MPIO_Finalize_async_io_thread: end\n");

#endif /* MPICH_THREAD_LEVEL == MPI_THREAD_MULTIPLE */
    return mpi_errno;
}
