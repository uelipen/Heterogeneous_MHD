#ifndef H_DMA_LIST_ELEM_H 
#define H_DMA_LIST_ELEM_H

struct dma_list_elem {
	union {
		unsigned int all32;
		struct {
			unsigned int stall: 	1;
			unsigned int reserved:	15;
			unsigned int nbytes:	16;
		} bits;
	} size;
	unsigned int ea_low;
};


#endif

