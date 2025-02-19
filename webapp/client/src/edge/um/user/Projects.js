import React from 'react'
import { Button } from 'reactstrap'
import { useNavigate } from 'react-router-dom'
import ProjectTable from '../common/ProjectTable'

const Projects = (props) => {
  const navigate = useNavigate()
  return (
    <div className="animated fadeIn">
      <div className="edge-right">
        <Button
          color="primary"
          size="sm"
          className="rounded-pill"
          onClick={() => navigate('/user/allProjects')}
          outline
        >
          &nbsp;Show all projects available to me&nbsp;
        </Button>
      </div>
      <br></br>
      <ProjectTable tableType="user" title={'My Projects'} {...props} />
      <br></br>
    </div>
  )
}

export default Projects
