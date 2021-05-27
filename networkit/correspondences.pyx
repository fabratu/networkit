# distutils: language=c++

from libc.stdint cimport uint64_t

ctypedef uint64_t count
ctypedef uint64_t index

from .structures cimport _Partition, Partition

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

cdef extern from "<networkit/correspondences/Correspondences.hpp>":
	cdef cppclass _Correspondences "NetworKit::Correspondences":
		_Correspondences() except +
		count run(_Partition, _Partition) except +

cdef class Correspondences:
	""" Finds correspondences between two partitions and provides details on correspondences
	"""
	cdef _Correspondences* _this

	def __cinit__(self):
		self._this = new _Correspondences()

	def run(self, Partition partition1, Partition partition2):
		return self._this.run(partition1._this, partition2._this)
