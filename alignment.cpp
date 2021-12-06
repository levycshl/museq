#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#define EP(A,i,j)   (*(double *)((A)->data + (i)*(A)->strides[0] + (j)*(A)->strides[1]))
#define EP2I(A,i,j) (*(long int *)   ((A)->data + (i)*(A)->strides[0] + (j)*(A)->strides[1]))
#define EP1I(A,i)   (*(long int *)   ((A)->data + (i)*(A)->strides[0]))


#define LEFT 0
#define DIAG 1
#define UP   2
#define INITIALIZE 0


static PyObject *
needleman(PyObject *self, PyObject *args)
{
    // inputs are:
    // S1, S2       integer [N1], [N2]
    // score_matrix integer [4x4]            
    // gap_costs    integer [2]
	//
	// empty return values:
    // nw_score     integer [N1 + 1, N2 + 1]
    // nw_trace     integer [N1 + 1, N2 + 1]
    const PyArrayObject * S1;
    const PyArrayObject * S2;
    const PyArrayObject * score_matrix;
    const PyArrayObject * gap_costs;
    const PyArrayObject *nw_score;
    const PyArrayObject *nw_trace;


	int score_left, score_diag, score_up;

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!", 
			&PyArray_Type, &S1, 
			&PyArray_Type, &S2,
			&PyArray_Type, &score_matrix, 
			&PyArray_Type, &gap_costs, 
			&PyArray_Type, &nw_score, 
			&PyArray_Type, &nw_trace))
       return NULL;

    const long int N1{S1->dimensions[0]};
    const long int N2{S2->dimensions[0]};
    const long int gap_start{EP1I(gap_costs, 0)};
    const long int gap_extend{EP1I(gap_costs, 1)};
	// initialize the score matrix with the gaps
	if(INITIALIZE){
	    int present_value = gap_start;
	    for (int i=1; i<=N1; i++) {
			EP2I(nw_score, i, 0) = present_value;
			present_value += gap_extend;
	    	}
		present_value = gap_start;
		for (int i=1; i<=N2; i++) {
			EP2I(nw_score, 0, i) = present_value;
			present_value += gap_extend;
			}
	}
	// fill in the matricies
	for (int i=0; i< N1; i++) {
		for(int j=0; j<N2; j++) {
			score_left = EP2I(nw_score, i+1, j  );
			score_up   = EP2I(nw_score, i  , j+1);
			score_diag = EP2I(nw_score, i  , j  );
			// adjust penalties depending on affine gaps and match_score
			score_left += ( EP2I(nw_trace, i+1, j  ) == LEFT ) ? gap_extend : gap_start;
			score_up   += ( EP2I(nw_trace, i  , j+1) == UP   ) ? gap_extend : gap_start;
			score_diag += EP2I(score_matrix, EP1I(S1, i), EP1I(S2, j));
			if (score_diag >= score_left){
				if (score_diag >= score_up) {
					EP2I(nw_trace, i+1, j+1) = DIAG;
					EP2I(nw_score, i+1, j+1) = score_diag;
					}
				else {
					EP2I(nw_trace, i+1, j+1) = UP;
					EP2I(nw_score, i+1, j+1) = score_up;
					}
			}
			else {
				if(score_left >= score_up){
					EP2I(nw_trace, i+1, j+1) = LEFT;
					EP2I(nw_score, i+1, j+1) = score_left;
					}
				else {
					EP2I(nw_trace, i+1, j+1) = UP;
					EP2I(nw_score, i+1, j+1) = score_up;
					}
			}					
		}
	}	
    return Py_BuildValue("i", sizeof(long int));
}

static PyMethodDef AlignmentMethods[] = {
    {"needleman",  needleman, METH_VARARGS, "Perform needleman-wunsch"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initAlignmentExt(void)
{
    (void) Py_InitModule("AlignmentExt", AlignmentMethods);
    import_array(); 
}

