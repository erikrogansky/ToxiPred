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

# Production commands
prod-up:
	docker compose -f docker-compose.prod.yml up -d

prod-down:
	docker compose -f docker-compose.prod.yml down

prod-logs:
	docker compose -f docker-compose.prod.yml logs -f

prod-status:
	./scripts/deploy.sh status

deploy:
	./scripts/deploy.sh deploy

rollback:
	./scripts/deploy.sh rollback

backup:
	./scripts/migrate.sh

# Utility commands
clean:
	docker system prune -f

clean-all:
	docker system prune -af
	docker volume prune -f

# Health & monitoring
health:
	./scripts/health-check.sh
