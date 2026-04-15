/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

//#include "mpiimpl.h"
//#include "mpiu_thread.h"

//#include "mpioimpl.h"
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "mpioimpl.h"
#include "mpiu_greq.h"
#include "async_io_thread.h"
#include "async_io_funcs.h"

/*
 *   Global varibles
 */

// ToDo: Move to FLexMPI library
long int    __attribute__ ((visibility ("default"))) EMPI_IOBLOCK = 0;
long int    __attribute__ ((visibility ("default"))) EMPI_DATA_READ = 0;
long int    __attribute__ ((visibility ("default"))) EMPI_DATA_IREAD = 0;
long int    __attribute__ ((visibility ("default"))) EMPI_DATA_WRITE = 0;
long int    __attribute__ ((visibility ("default"))) EMPI_DATA_IWRITE = 0;
long int    __attribute__ ((visibility ("default"))) EMPI_UTIME_READ = 0;
long int    __attribute__ ((visibility ("default"))) EMPI_UTIME_IREAD = 0;
long int    __attribute__ ((visibility ("default"))) EMPI_UTIME_WRITE = 0;
long int    __attribute__ ((visibility ("default"))) EMPI_UTIME_IWRITE = 0;
long double __attribute__ ((visibility ("default"))) EMPI_DESIRED_BW_READ = 0.0;
long double __attribute__ ((visibility ("default"))) EMPI_DESIRED_BW_IREAD = 0.0;
long double __attribute__ ((visibility ("default"))) EMPI_DESIRED_BW_WRITE = 0.0;
long double __attribute__ ((visibility ("default"))) EMPI_DESIRED_BW_IWRITE = 0.0;
long double __attribute__ ((visibility ("default"))) EMPI_SCALE_BW_READ = 1.0;
long double __attribute__ ((visibility ("default"))) EMPI_SCALE_BW_IREAD = 1.0;
long double __attribute__ ((visibility ("default"))) EMPI_SCALE_BW_WRITE = 1.0;
long double __attribute__ ((visibility ("default"))) EMPI_SCALE_BW_IWRITE = 1.0;
int         __attribute__ ((visibility ("default"))) EMPI_WORLD_RANK = 0;
int         __attribute__ ((visibility ("default"))) EMPI_WORLD_SIZE = 1;

/*
 *   local functions
 */

/* Starts a generalized request  */

static void generic_request_create_start(MPI_File *fh, MPI_Offset bytes,
        int *error_code, MPI_Request *request)
{
    MPI_Status *status;
    status = (MPI_Status *)ADIOI_Malloc(sizeof(MPI_Status));

    status->MPI_ERROR = *error_code;
#ifdef HAVE_STATUS_SET_BYTES
    MPIR_Status_set_bytes(status, MPI_BYTE, bytes);
#endif
    /* --BEGIN ERROR HANDLING-- */
    if (*error_code != MPI_SUCCESS)
        *error_code = MPIO_Err_return_file(*fh, *error_code);
    /* --END ERROR HANDLING-- */
    MPI_Grequest_start(MPIU_Greq_query_fn, MPIU_Greq_free_fn,
            MPIU_Greq_cancel_fn, status, request);
}

static void generic_request_complete(MPI_Request *request)
{
    MPI_Grequest_complete(*request);
}

/*@
    WriteContig - Wraper to ADIO_WriteContig
 @*/
static void WriteContig (ADIO_File fd, const void *buf, int count, MPI_Datatype datatype, int file_ptr_type, ADIO_Offset offset, useconds_t exec_time, MPI_Status *status, int *error_code, int *time_excess)
{
    double ltime_ini;
    double ltime_middle;
    double ltime_end;
    int time_wait;
    ltime_ini = MPI_Wtime ();
    ROMIO_THREAD_CS_ENTER();
    ADIO_WriteContig(fd, buf, count, datatype, file_ptr_type, offset, status, error_code);
    ROMIO_THREAD_CS_EXIT();
    ltime_middle = MPI_Wtime ();
    if (exec_time != 0) {
        int run_time = (int)((ltime_middle - ltime_ini)*1000000.0);
        time_wait = exec_time - run_time;
        PRINTF("WriteContig: exec_time=%d, run_time=%d, ltime_ini=%f, ltime_middle=%f\n",  exec_time, run_time, ltime_ini, ltime_middle);

        if ((time_wait < 0)) {
            (*time_excess) = (*time_excess) - time_wait;
            time_wait = 0;
        } else if ((time_wait > 0) & (time_wait < (*time_excess))) {
            (*time_excess) = (*time_excess) - time_wait;
            time_wait = 0;
        } else if (time_wait > (*time_excess)) {
            time_wait = time_wait - (*time_excess);
            (*time_excess) = 0;
            // Blocks a certain amount of time
            usleep(time_wait);
        }
    }
    ltime_end = MPI_Wtime ();
    int real_wait_time = (int)((ltime_end - ltime_middle)*1000000.0);
    (*time_excess) = (*time_excess) + (real_wait_time - time_wait);
}

/*@
    ReadContig - Wraper to ADIO_ReadContig
 @*/
static void ReadContig (ADIO_File fd, void *buf, int count, MPI_Datatype datatype, int file_ptr_type, ADIO_Offset offset, useconds_t exec_time, MPI_Status *status, int *error_code, int *time_excess)

{
    double ltime_ini;
    double ltime_middle;
    double ltime_end;
    int time_wait;

    ROMIO_THREAD_CS_ENTER();
    ltime_ini = MPI_Wtime ();
    ADIO_ReadContig(fd, buf, count, datatype, file_ptr_type, offset, status, error_code);
    ROMIO_THREAD_CS_EXIT();
    ltime_middle = MPI_Wtime ();
    if (exec_time != 0) {
        int run_time = (int)((ltime_middle - ltime_ini)*1000000.0);
        time_wait = exec_time - run_time;
        PRINTF("ReadContig: exec_time=%d, run_time=%d, ltime_ini=%f, ltime_middle=%f\n",  exec_time, run_time, ltime_ini, ltime_middle);

        if ((time_wait < 0)) {
            (*time_excess) = (*time_excess) - time_wait;
            time_wait = 0;
        } else if ((time_wait > 0) & (time_wait < (*time_excess))) {
            (*time_excess) = (*time_excess) - time_wait;
            time_wait = 0;
        } else if (time_wait > (*time_excess)) {
            time_wait = time_wait - (*time_excess);
            (*time_excess) = 0;
            // Blocks a certain amount of time
            usleep(time_wait);
        }
    }
    ltime_end = MPI_Wtime ();
    int real_wait_time = (int)((ltime_end - ltime_middle)*1000000.0);
    (*time_excess) = (*time_excess) + (real_wait_time - time_wait);
}

/*@
    IwriteContig - Wraper to ADIO_WriteContig
 @*/
static void IwriteContig (ADIO_File fd, const void *buf, int count, MPI_Datatype datatype, int file_ptr_type, ADIO_Offset offset, useconds_t exec_time, MPI_Status *status, int *error_code, int *time_excess)


{
    WriteContig (fd, buf, count, datatype, file_ptr_type, offset, exec_time, status, error_code, time_excess);
}

/*@
    IreadContig - Wraper to ADIO_ReadContig
 @*/
static void IreadContig (ADIO_File fd, void *buf, int count, MPI_Datatype datatype, int file_ptr_type, ADIO_Offset offset, useconds_t exec_time, MPI_Status *status, int *error_code, int *time_excess)


{
    ReadContig (fd, buf, count, datatype, file_ptr_type, offset, exec_time, status, error_code, time_excess);
}

/*@
 exec_operation - select file operaation for bw control
 @*/
static void exec_operation(void *buf, int count, MPI_Datatype datatype, ADIO_Offset offset, int exec_time, int oper, void *oper_data, int oper_data_size, MPI_Status *status, int *error_code, int *time_excess)
{
    PRINTF ("exec_operation: begin\n");
    if (oper == ASYNC_OP_WriteContig) {
        Async_WriteContig_data_t * cast_data;
        cast_data = (Async_WriteContig_data_t *)oper_data;
        PRINTF("exec_operation: ASYNC_OP_WriteContig: exec_time=%d, time_excess=%d, fd=%p, buf/orig_buf=%p/%p count/orig_count=%d/%d, dataype/orig_dataype=%d/%d, file_ptr_type=%d, offset/orig_offset=%lld/%lld, status/orig_status=%p/%p, error_code/orig_error_code=%p/%p\n", exec_time, (*time_excess), cast_data->fd, buf, cast_data->buf, count, cast_data->count, datatype,  cast_data->datatype, cast_data->file_ptr_type, offset, cast_data->offset, status, cast_data->status, error_code, cast_data->error_code);
        WriteContig(cast_data->fd, buf, count, datatype, cast_data->file_ptr_type, offset, exec_time, status, error_code, time_excess);

    } else if (oper == ASYNC_OP_ReadContig) {
        Async_ReadContig_data_t * cast_data;
        cast_data = (Async_ReadContig_data_t *)oper_data;
        PRINTF("exec_operation: ASYNC_OP_ReadContig: exec_time=%d, time_excess=%d, fd=%p, buf/orig_buf=%p/%p count/orig_count=%d/%d, dataype/orig_dataype=%d/%d, file_ptr_type=%d, offset/orig_offset=%lld/%lld, status/orig_status=%p/%p, error_code/orig_error_code=%p/%p\n", exec_time, (*time_excess), cast_data->fd, buf, cast_data->buf, count, cast_data->count, datatype,  cast_data->datatype, cast_data->file_ptr_type, offset, cast_data->offset, status, cast_data->status, error_code, cast_data->error_code);
        ReadContig(cast_data->fd, buf, count, datatype, cast_data->file_ptr_type, offset, exec_time, status, error_code, time_excess);

    } else if (oper == ASYNC_OP_IWriteContig) {
        Async_IWriteContig_data_t * cast_data;
        cast_data = (Async_IWriteContig_data_t *)oper_data;
        PRINTF("exec_operation: ASYNC_OP_IWriteContig: exec_time=%d, time_excess=%d, fd=%p, buf/orig_buf=%p/%p count/orig_count=%d/%d, dataype/orig_dataype=%d/%d, file_ptr_type=%d, offset/orig_offset=%lld/%lld, status=%p, error_code=%p request=%p\n", exec_time, (*time_excess), cast_data->fd, buf, cast_data->buf, count, cast_data->count, datatype,  cast_data->datatype, cast_data->file_ptr_type, offset, cast_data->offset, status, error_code, cast_data->request);
        IwriteContig(cast_data->fd, buf, count, datatype, cast_data->file_ptr_type, offset, exec_time, status, error_code, time_excess);

    } else if (oper == ASYNC_OP_IReadContig) {
        Async_IReadContig_data_t * cast_data;
        cast_data = (Async_IReadContig_data_t *)oper_data;
        PRINTF("exec_operation: ASYNC_OP_IReadContig: exec_time=%d, time_excess=%d, fd=%p, buf/orig_buf=%p/%p count/orig_count=%d/%d, dataype/orig_dataype=%d/%d, file_ptr_type=%d, offset/orig_offset=%lld/%lld, status=%p, error_code=%p request=%p\n", exec_time, (*time_excess), cast_data->fd, buf, cast_data->buf, count, cast_data->count, datatype,  cast_data->datatype, cast_data->file_ptr_type, offset, cast_data->offset, status, error_code, cast_data->request);
        IreadContig(cast_data->fd, buf, count, datatype, cast_data->file_ptr_type, offset, exec_time, status, error_code, time_excess);

    } else {
        PRINTF ("exec_operation: ERROR oper=%d not valid\n",oper);
    }
    PRINTF ("exec_operation: end\n");
}
/*@
 control_bw_operation - implement generic bw control
 @*/
static void control_bw_operation (void *buf, int count, MPI_Datatype datatype, ADIO_Offset offset, int oper, void *oper_data, int oper_data_size)
{
    int async_error_code;
    int time_excess=0;
    double time_ini=0;
    double time_end=0;
    long int aggr_data=0;
    long double* EMPI_SCALE_BW = NULL;
    long double* EMPI_DESIRED_BW = NULL;
    int datatype_size;

    PRINTF("control_bw_operation: begin\n");

    // get initial time for calculate the bandwitdh
    time_ini = MPI_Wtime ();

    // get the size of the datatype
    MPI_Type_size(datatype,&datatype_size);

    if (oper == ASYNC_OP_ReadContig) {
        EMPI_SCALE_BW = &EMPI_SCALE_BW_READ;
        EMPI_DESIRED_BW = &EMPI_DESIRED_BW_READ;
    } else if (oper == ASYNC_OP_IReadContig) {
        EMPI_SCALE_BW = &EMPI_SCALE_BW_IREAD;
        EMPI_DESIRED_BW = &EMPI_DESIRED_BW_IREAD;
    } else if (oper == ASYNC_OP_WriteContig) {
        EMPI_SCALE_BW = &EMPI_SCALE_BW_WRITE;
        EMPI_DESIRED_BW = &EMPI_DESIRED_BW_WRITE;
    } else if (oper == ASYNC_OP_IWriteContig) {
        EMPI_SCALE_BW = &EMPI_SCALE_BW_IWRITE;
        EMPI_DESIRED_BW = &EMPI_DESIRED_BW_IWRITE;
    }
    
    PRINTF("control_bw_operation: --> EMPI_WORLD_RANK: %d, EMPI_WORLD_SIZE: %d, EMPI_IOBLOCK: %ld, EMPI_DATA_READ: %ld, EMPI_UTIME_READ: %ld, EMPI_DATA_IREAD: %ld, EMPI_UTIME_IREAD: %ld, EMPI_DATA_WRITE: %ld, EMPI_UTIME_WRITE: %ld, EMPI_DATA_IWRITE: %ld, EMPI_UTIME_IWRITE: %ld, EMPI_DESIRED_BW: %f, EMPI_SCALE_BW: %f\n",EMPI_WORLD_RANK, EMPI_WORLD_SIZE, EMPI_IOBLOCK, EMPI_DATA_READ, EMPI_UTIME_READ, EMPI_DATA_IREAD, EMPI_UTIME_IREAD, EMPI_DATA_WRITE,  EMPI_UTIME_WRITE, EMPI_DATA_IWRITE,  EMPI_UTIME_IWRITE, EMPI_DESIRED_BW, EMPI_SCALE_BW);
    
    if(*EMPI_DESIRED_BW!=0.0){
        // IO Scheduling variables
        double adjusted_desired_bw=0;
        
        int exec_time_block = 0;
        int exec_time = 0;

        int num_phases,i;
        long int array_offset,local_count;
        MPI_Offset local_off;
        
        // ToDo: include non-multiple of I/O phases
        num_phases=(int)ceil((double)count/(double)EMPI_IOBLOCK);
        
        // Get adjusted_desired_bw
        

        // Get exec time for a standard block
        //exec_time_block = (int)((((double)(EMPI_WORLD_SIZE*EMPI_IOBLOCK*datatype_size)) / ((adjusted_desired_bw)*1024.0*1024.0))*1000000.0);
        //exec_time_block = (int)((((double)(EMPI_IOBLOCK*datatype_size)) / ((adjusted_desired_bw)*1024.0*1024.0))*1000000.0);
        
        PRINTF("control_bw_operation: exec_time_block: %d, EMPI_DESIRED_BW: %lf, adjusted_desired_bw: %lf, EMPI_WORLD_SIZE: %d, EMPI_IOBLOCK: %ld, datatype_size:%d\n",exec_time_block,EMPI_DESIRED_BW, adjusted_desired_bw, EMPI_WORLD_SIZE, EMPI_IOBLOCK, datatype_size);

        for (i=0;i<num_phases;i++){

            adjusted_desired_bw = (double)(*EMPI_DESIRED_BW) * (double)(*EMPI_SCALE_BW);
            exec_time_block = (int)((((double)(EMPI_IOBLOCK*datatype_size)) / (adjusted_desired_bw))*1000000.0);


            local_off=offset+(MPI_Offset)(i*((MPI_Offset)(EMPI_IOBLOCK*(long int)datatype_size)));
            array_offset=((long int)i)*EMPI_IOBLOCK;
            if(i<num_phases-1){
                local_count=EMPI_IOBLOCK;
                exec_time = exec_time_block;
            }
            else
            {
                local_count=((long int)count)-EMPI_IOBLOCK*(long int)i;
                //exec_time = (int)((((double)(EMPI_WORLD_SIZE*local_count*datatype_size)) / ((adjusted_desired_bw)*1024.0*1024.0))*1000000.0);
                //exec_time = (int)((((double)(local_count*datatype_size)) / ((adjusted_desired_bw)*1024.0*1024.0))*1000000.0);
                exec_time = (int)((((double)(local_count*datatype_size)) / (adjusted_desired_bw))*1000000.0);
            }
            PRINTF("control_bw_operation:  ----> Rank [%d] i: %d count: %ld  offset: %ld  array_offset: %ld --- Block: %ld - exec_time: %d \n",EMPI_WORLD_RANK,i,(long int)local_count,(long int)local_off,(long int)array_offset,EMPI_IOBLOCK,exec_time);
            
            exec_operation((void *)buf+(int)array_offset*datatype_size, (int)local_count, datatype, local_off, (useconds_t)exec_time, oper, oper_data, oper_data_size, MPI_STATUS_IGNORE, &async_error_code, &time_excess);

            aggr_data+=(local_count*datatype_size);
            PRINTF ("control_bw_operation: --> aggr_data: %ld, local_count: %ld, datatype_size: %d,\n",aggr_data, local_count, datatype_size);
            
         }
                 
    }
    else{
        PRINTF("control_bw_operation: --> EMPI_WORLD_RANK: %d, EMPI_WORLD_SIZE: %d, EMPI_IOBLOCK: %ld, EMPI_DATA_READ: %ld, EMPI_UTIME_READ: %ld, EMPI_DATA_IREAD: %ld, EMPI_UTIME_IREAD: %ld, EMPI_DATA_WRITE: %ld, EMPI_UTIME_WRITE: %ld, EMPI_DATA_IWRITE: %ld, EMPI_UTIME_IWRITE: %ld, EMPI_DESIRED_BW: %f, EMPI_SCALE_BW: %f\n",EMPI_WORLD_RANK, EMPI_WORLD_SIZE, EMPI_IOBLOCK, EMPI_DATA_READ, EMPI_UTIME_READ, EMPI_DATA_IREAD, EMPI_UTIME_IREAD, EMPI_DATA_WRITE,  EMPI_UTIME_WRITE, EMPI_DATA_IWRITE,  EMPI_UTIME_IWRITE, EMPI_DESIRED_BW, EMPI_SCALE_BW);
        exec_operation(buf, count, datatype, offset, 0, oper, oper_data, oper_data_size, MPI_STATUS_IGNORE, &async_error_code, &time_excess);
        aggr_data=(count*datatype_size);

    }

    // get initial time for calculate the bandwitdh
    time_end = MPI_Wtime ();
    
    if (oper == ASYNC_OP_ReadContig) {
        EMPI_SCALE_BW_READ = (long double)*EMPI_SCALE_BW;
        EMPI_DATA_READ = EMPI_DATA_READ + aggr_data;
        EMPI_UTIME_READ = EMPI_UTIME_READ + (long int)((time_end - time_ini) * 1000000.0) ;
    } else if  (oper == ASYNC_OP_IReadContig) {
        EMPI_SCALE_BW_IREAD = (long double)*EMPI_SCALE_BW;
        EMPI_DATA_IREAD = EMPI_DATA_IREAD + aggr_data;
        EMPI_UTIME_IREAD = EMPI_UTIME_IREAD + (long int)((time_end - time_ini) * 1000000.0) ;
    } else if (oper == ASYNC_OP_WriteContig) {
        EMPI_SCALE_BW_WRITE = (long double)*EMPI_SCALE_BW;
        EMPI_DATA_WRITE = EMPI_DATA_WRITE + aggr_data;
        EMPI_UTIME_WRITE = EMPI_UTIME_WRITE + (long int)((time_end - time_ini) * 1000000.0) ;
    } else if (oper == ASYNC_OP_IWriteContig) {
        EMPI_SCALE_BW_IWRITE = (long double)*EMPI_SCALE_BW;
        EMPI_DATA_IWRITE = EMPI_DATA_IWRITE + aggr_data;
        EMPI_UTIME_IWRITE = EMPI_UTIME_IWRITE + (long int)((time_end - time_ini) * 1000000.0) ;
    }
    PRINTF("control_bw_operation:  --> Rank %d    EMPI_DATA_READ: %ld, EMPI_UTIME_READ: %ld, EMPI_DATA_IREAD: %ld, EMPI_UTIME_IREAD: %ld, EMPI_DATA_WRITE: %ld, EMPI_UTIME_WRITE: %ld, EMPI_DATA_IWRITE: %ld, EMPI_UTIME_IWRITE: %ld, time_ini=%lf time_end=%lf\n ",EMPI_WORLD_RANK, EMPI_DATA_READ, EMPI_UTIME_READ, EMPI_DATA_IREAD, EMPI_UTIME_IREAD, EMPI_DATA_WRITE,  EMPI_UTIME_WRITE, EMPI_DATA_IWRITE,  EMPI_UTIME_IWRITE, time_ini, time_end);

    PRINTF("control_bw_operation: end\n");
}


/*@
    Async_exec_op - exec a thread-based async op

Input Parameters:
 . oper - operation to execute
 . data - data for the operation
 . size - size of the data buffer

.N fortran modificar
@*/

int Async_exec_op(int oper, void *data, int size)
{
    PRINTF ("Async_exec_op: begin\n");
    if (oper == ASYNC_OP_WriteContig) {
        Async_WriteContig_data_t * cast_data;
        cast_data = (Async_WriteContig_data_t *)data;
        PRINTF("Async_exec_op: ASYNC_OP_WriteContig: fd=%p, buf=%p, count=%d, dataype=%d, file_ptr_type=%d, offset=%lld, status=%p, error_code=%p\n", cast_data->fd, cast_data->buf, cast_data->count, cast_data->datatype, cast_data->file_ptr_type, cast_data->offset, cast_data->status, cast_data->error_code);
        
        control_bw_operation ((void *)(cast_data->buf), cast_data->count, cast_data->datatype, cast_data->offset, oper, data, size);

    } else if (oper == ASYNC_OP_ReadContig) {
            Async_ReadContig_data_t * cast_data;
            cast_data = (Async_ReadContig_data_t *)data;
            PRINTF("Async_exec_op: ASYNC_OP_ReadContig: fd=%p, buf=%p, count=%d, dataype=%d, file_ptr_type=%d, offset=%lld, status=%p, error_code=%p\n", cast_data->fd, cast_data->buf, cast_data->count, cast_data->datatype, cast_data->file_ptr_type, cast_data->offset, cast_data->status, cast_data->error_code);
            
            control_bw_operation (cast_data->buf, cast_data->count, cast_data->datatype, cast_data->offset, oper, data, size);

    } else if (oper == ASYNC_OP_IWriteContig) {
        Async_IWriteContig_data_t * cast_data;
        cast_data = (Async_IWriteContig_data_t *)data;
        PRINTF("Async_exec_op: ASYNC_OP_IWriteContig: fd=%p, buf=%p, count=%d, dataype=%d, file_ptr_type=%d, offset=%lld, request=%p\n", cast_data->fd, cast_data->buf, cast_data->count, cast_data->datatype, cast_data->file_ptr_type, cast_data->offset, cast_data->request);

        control_bw_operation ((void *)(cast_data->buf), cast_data->count, cast_data->datatype, cast_data->offset, oper, data, size);

        // signal the end of the write operation
        generic_request_complete(cast_data->request);
    
    } else if (oper == ASYNC_OP_IReadContig) {
        Async_IReadContig_data_t * cast_data;
        cast_data = (Async_IReadContig_data_t *)data;
        PRINTF("Async_exec_op: ASYNC_OP_IReadContig: fd=%p, buf=%p, count=%d, dataype=%d, file_ptr_type=%d, offset=%lld, request=%p\n", cast_data->fd, cast_data->buf, cast_data->count, cast_data->datatype, cast_data->file_ptr_type, cast_data->offset, cast_data->request);

        control_bw_operation (cast_data->buf, cast_data->count, cast_data->datatype, cast_data->offset, oper, data, size);

        // signal the end of the read operation
        generic_request_complete(cast_data->request);

    } else if (oper == ASYNC_OP_IWriteContig) {
        PRINTF ("Async_exec_op: Wait oper\n");
    } else {
        PRINTF ("Async_exec_op: ERROR oper=%d not valid\n",oper);
    }
    PRINTF ("Async_exec_op: end\n");
    return 0;
}

/*@
    Async_WriteContig - use a thread-based ADIO_WriteContig
 @*/

void Async_WriteContig(ADIO_File fd, const void *buf, int count, MPI_Datatype datatype, int file_ptr_type, ADIO_Offset offset, MPI_Status *status, int *error_code)
{
    PRINTF ("Async_WriteContig: begin\n");
    ROMIO_THREAD_CS_EXIT();
    Async_WriteContig_data_t iwritecontig_data = {fd,buf,count,datatype,file_ptr_type,offset,status,error_code};
    MPIIO_Enqueue_async_io_queue(&MPIIO_Async_io_queue, ASYNC_OP_WriteContig, (char *)&iwritecontig_data, sizeof(Async_WriteContig_data_t));
    Async_dummy_data_t wait_data = {0};
    MPIIO_Enqueue_async_io_queue(&MPIIO_Async_io_queue, ASYNC_OP_WAIT, (char *)&wait_data, sizeof(Async_dummy_data_t));
    MPIIO_WaitEmpty_async_io_queue(&MPIIO_Async_io_queue);
    ROMIO_THREAD_CS_ENTER();
    PRINTF ("Async_WriteContig: end\n");
}
 
/*@
    Async_ReadContig - use a thread-based ADIO_ReadContig
 @*/

void Async_ReadContig(ADIO_File fd, void *buf, int count, MPI_Datatype datatype, int file_ptr_type, ADIO_Offset offset, MPI_Status *status, int *error_code)
{
    PRINTF ("Async_ReadContig: begin\n");
    ROMIO_THREAD_CS_EXIT();
    Async_ReadContig_data_t ireadcontig_data = {fd,buf,count,datatype,file_ptr_type,offset,status,error_code};
    MPIIO_Enqueue_async_io_queue(&MPIIO_Async_io_queue, ASYNC_OP_ReadContig, (char *)&ireadcontig_data, sizeof(Async_ReadContig_data_t));
    Async_dummy_data_t wait_data = {0};
    MPIIO_Enqueue_async_io_queue(&MPIIO_Async_io_queue, ASYNC_OP_WAIT, (char *)&wait_data, sizeof(Async_dummy_data_t));
    MPIIO_WaitEmpty_async_io_queue(&MPIIO_Async_io_queue);
    ROMIO_THREAD_CS_ENTER();
    PRINTF ("Async_ReadContig: end\n");
}

/*@
    Async_IwriteContig - use a thread-based ADIO_IwriteContig
 @*/

void Async_IwriteContig(ADIO_File fd, const void *buf, int count, MPI_Datatype datatype, int file_ptr_type, ADIO_Offset offset, ADIO_Request *request, int *error_code)
{
    PRINTF ("Async_IwriteContig: begin\n");
    MPI_Count datatype_size;
    MPI_Type_size_x(datatype, &datatype_size);
    MPI_Offset nbytes = count * datatype_size;

    (*error_code) = MPI_SUCCESS;
    generic_request_create_start(&fd, nbytes, error_code, request);
    Async_IWriteContig_data_t elem = {fd,buf,count,datatype,file_ptr_type,offset,request};
    ROMIO_THREAD_CS_EXIT();
    MPIIO_Enqueue_async_io_queue(&MPIIO_Async_io_queue, ASYNC_OP_IWriteContig, (char *)&elem, sizeof(Async_IWriteContig_data_t));
    ROMIO_THREAD_CS_ENTER();
    //MPIIO_WaitEmpty_async_io_queue(&MPIIO_Async_io_queue);

    PRINTF ("Async_IwriteContig: end\n");
}

/*@
    Async_IreadContig - use a thread-based ADIO_IreadContig
 @*/

void Async_IreadContig(ADIO_File fd, void *buf, int count, MPI_Datatype datatype, int file_ptr_type, ADIO_Offset offset, ADIO_Request *request, int *error_code)
{
    PRINTF ("Async_IreadContig: begin\n");
    MPI_Count datatype_size;
    MPI_Type_size_x(datatype, &datatype_size);
    MPI_Offset nbytes = count * datatype_size;

    (*error_code) = MPI_SUCCESS;
    generic_request_create_start(&fd, nbytes, error_code, request);
    Async_IReadContig_data_t elem = {fd,buf,count,datatype,file_ptr_type,offset,request};
    MPIIO_Enqueue_async_io_queue(&MPIIO_Async_io_queue, ASYNC_OP_IReadContig, (char *)&elem, sizeof(Async_IReadContig_data_t));
    //MPIIO_WaitEmpty_async_io_queue(&MPIIO_Async_io_queue);

    PRINTF ("Async_IreadContig: end\n");
}

