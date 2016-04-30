#ifndef SPU_MFCIO_EXT_H 
#define SPU_MFCIO_EXT_H 

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdint.h>
#include <stdio.h>
#include <ctype.h>

static uint32_t msg[4]__attribute__ ((aligned (16)));

//	signal notify 1 and 2 registers definitions
#define SPU_SIG_NOTIFY_OFFSET 0x0C // offset from signal areas base
#define SPU_SIG_NOTIFY_OFFSET_SLOT 0x3 // 16B alignment (OFFSET&0xF)>>2 

// signal a remote SPE<92>s signal1 register
inline int write_signal1(uint32_t data, uint64_t ea_sig1, uint32_t tag_id)
{
        uint64_t ea_sig1_notify = ea_sig1 + SPU_SIG_NOTIFY_OFFSET;
        uint32_t idx;

        //printf("<SPE: write_signal1: starts\n" );

        //printf("<SPE: write_signal1:  ea_mfc=0x%llx, ea_in_mbox=0x%llx\n", ea_sig1, ea_sig1_notify );

        idx = SPU_SIG_NOTIFY_OFFSET_SLOT;
        msg[idx] = data;

        mfc_sndsig( &msg[idx], ea_sig1_notify, tag_id, 0,0)     ;
        mfc_write_tag_mask(1<<tag_id);
//printf("in signal1\n");
        mfc_read_tag_status_all();
        //mfc_read_tag_status_any();

        //printf("<SPE: write_in_mbox: complete\n" );

        return 1; // number of mailbox being written
}

// signal a remote SPE<92>s signal1 register
inline int write_signal2(uint32_t data, uint64_t ea_sig2, uint32_t tag_id)
{
        uint64_t ea_sig2_notify = ea_sig2 + SPU_SIG_NOTIFY_OFFSET;
        uint32_t idx;

        //printf("<SPE: write_signal2: starts\n" );

        //printf("<SPE: write_signal2:  ea_mfc=0x%llx, ea_in_mbox=0x%llx\n", ea_sig2, ea_sig2_notify );

        idx = SPU_SIG_NOTIFY_OFFSET_SLOT;
        msg[idx] = data;

        mfc_sndsig( &msg[idx], ea_sig2_notify, tag_id, 0,0)     ;
        mfc_write_tag_mask(1<<tag_id);
        mfc_read_tag_status_any();

        //printf("<SPE: write_in_mbox: complete\n" );

        return 1; // number of mailbox being written
}

#endif
