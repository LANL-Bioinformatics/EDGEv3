# Development

## Dependencies

* node.js v20
* MongoDB

## Install the webapp

    cd installation
    ./install.sh

## Configure and start ui client

    cd webapp/client
    cp .env.example .env
    (update settings in .env)
    npm start
    
## Configure and start api server

    cd webapp/server
    cp .env.example .env
    (update settings in .env)
    npm start


## View the website

    http://localhost:3000

## Note

- Have to restart the api server when any changes made in webapp/server code or webapp/server/.env.
    