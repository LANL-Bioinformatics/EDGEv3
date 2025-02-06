import React from 'react'
import { Badge } from 'reactstrap'
import ProjectTable from '../common/ProjectTable'

const Projects = (props) => {
  return (
    <div className="animated fadeIn">
      <div className="clearfix">
        <Badge color="danger" pill>
          Admin tool
        </Badge>
      </div>
      <br></br>
      <ProjectTable tableType="admin" title={'Manage Projects'} {...props} />
      <br></br>
    </div>
  )
}

export default Projects
