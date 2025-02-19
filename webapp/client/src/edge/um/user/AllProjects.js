import React from 'react'
import ProjectTableViewOnly from '../common/ProjectTableViewOnly'
import { apis } from '../../common/util'

const AllProjects = (props) => {
  return (
    <div className="animated fadeIn">
      <ProjectTableViewOnly
        title={'All Projects Available to Me'}
        projectPageUrl={'/user/project?code='}
        api={apis.userAllProjects}
      />
      <br></br>
    </div>
  )
}

export default AllProjects
