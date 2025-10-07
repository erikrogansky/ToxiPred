.PHONY: up down logs rebuild

build:
\tdocker compose up --build

up:
\tdocker compose up

down:
\tdocker compose down

logs:
\tdocker compose logs -f

server-logs:
\tdocker compose logs server-toxipred

app-logs:
\tdocker compose logs app-toxipred

rebuild:
\tdocker compose build --no-cache
