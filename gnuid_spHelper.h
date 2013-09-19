#include "libmesh/libmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparsity_pattern.h"

class AugmentSparsityPatternDefectiveBC : public DofMap::AugmentSparsityPattern
{
public:
  AugmentSparsityPatternDefectiveBC (MeshBase& mesh, DofMap& dofMap, EquationSystems& es)
  {
    _mesh = &mesh;
    _dofMap = & dofMap;
    _es = & es;
  }
  
private:
  void augment_sparsity_pattern (SparsityPattern::Graph & sparsity,
                                  std::vector<dof_id_type> & n_nz,
                                  std::vector<dof_id_type> & n_oz) 
  {
    unsigned int n_flow_bc = _es->parameters.get<unsigned int>("n flow bound");
    std::vector<std::vector<const Elem*> > bcElems(n_flow_bc);
    MeshBase::const_element_iterator el = _mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _mesh->active_local_elements_end(); 
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      for (unsigned int side=0; side<elem->n_sides(); side++)
      {
        if (elem->neighbor(side) == NULL)
        {
          const unsigned short boundary_id = _mesh->boundary_info->boundary_id(elem,side);
          for (unsigned int lid = 0; lid < n_flow_bc; lid++)
          {
            char flow_bc_id[1024];
            sprintf(flow_bc_id,"flow_bc_%02d_id",lid);
            if (boundary_id == _es->parameters.get<unsigned short>(flow_bc_id))
            {
	      bcElems[lid].push_back(elem);
	    }
	  }
	}
      }
    }

    const processor_id_type proc_id           = _mesh->processor_id();
    const dof_id_type n_dofs_on_proc    = _dofMap->n_dofs_on_processor(proc_id);
    const dof_id_type first_dof_on_proc = _dofMap->first_dof(proc_id);
    const dof_id_type end_dof_on_proc   = _dofMap->end_dof(proc_id);

    std::vector<dof_id_type> elementDofs, coupledElementDofs; 

    for (unsigned int lid = 0; lid < n_flow_bc; lid++)
      for (unsigned int i = 0; i < bcElems[lid].size(); i++)
      {
        const Elem* elem = bcElems[lid][i];
        _dofMap->dof_indices (elem, elementDofs);
        const unsigned int n_dofs_on_element = libmesh_cast_int<unsigned int>(elementDofs.size());

        for (unsigned int i=0; i<n_dofs_on_element; i++)
        {
          const dof_id_type ig = elementDofs[i];
          SparsityPattern::Row *row;
          
          if ((ig < first_dof_on_proc) || (ig >= end_dof_on_proc))
          {
            std::cout<<" error, dof not on proc";
          }
          row = &sparsity[ig - first_dof_on_proc];

          for (unsigned int j = 0; j < bcElems[lid].size(); j++)
          {
            const Elem* coupledElem = bcElems[lid][j];
            _dofMap->dof_indices (coupledElem, coupledElementDofs);

            const std::size_t n_dofs_on_coupled_element = coupledElementDofs.size();
            for (std::size_t j=0; j<n_dofs_on_coupled_element; j++)
            {
              const dof_id_type jg = coupledElementDofs[j];
              // See if jg is in the sorted range
              std::pair<SparsityPattern::Row::iterator,
                        SparsityPattern::Row::iterator>
                        pos = std::equal_range (row->begin(), row->end(), jg);
              // Insert jg if it wasn't found
              if (pos.first == pos.second)
                row->insert (pos.first, jg);
            }
          }
        }
      }
    for (dof_id_type i=0; i<n_dofs_on_proc; i++)
    {
      // Get the row of the sparsity pattern
      SparsityPattern::Row &row = sparsity[i];
      for (dof_id_type j=0; j<row.size(); j++)
        if ((row[j] < first_dof_on_proc) || (row[j] >= end_dof_on_proc))
          n_oz[i]++;
        else
          n_nz[i]++;
      row.clear();
    }
  }

  MeshBase * _mesh;
  DofMap * _dofMap;
  EquationSystems * _es;
};
