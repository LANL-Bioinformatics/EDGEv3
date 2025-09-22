FROM nginx:1.12

#  default conf for proxy service
COPY ./default.conf /etc/nginx/conf.d/default.conf
#COPY ./nginx.conf /etc/nginx/nginx.conf

# NOT FOUND response
COPY ./backend-not-found.html /var/www/html/backend-not-found.html

# Proxy configurations
COPY ./includes/ /etc/nginx/includes/
