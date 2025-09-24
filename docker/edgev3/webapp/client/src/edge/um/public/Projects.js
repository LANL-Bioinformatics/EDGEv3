import React from 'react'
import ProjectTableViewOnly from '../common/ProjectTableViewOnly'
import { apis } from '../../common/util'

const Projects = (props) => {
  return (
    <div className="animated fadeIn">
      <ProjectTableViewOnly
        title={'Public Projects'}
        projectPageUrl={'/public/project?code='}
        api={apis.publicProjects}
      />
      <br></br>
    </div>
  )
}

export default Projects
