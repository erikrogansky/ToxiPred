# frontend/Dockerfile
FROM node:20-alpine
WORKDIR /app

# bring the whole app so `quasar prepare` sees index.html, src/, etc.
COPY . .

RUN npm ci --no-audit --no-fund

EXPOSE 9000
ENV HOST=0.0.0.0
CMD ["npm", "run", "dev", "--", "--host", "0.0.0.0"]
