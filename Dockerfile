# ─────────────────────────────────────────────────────────────────────────────
# Stage 1: Build
# ─────────────────────────────────────────────────────────────────────────────
FROM mcr.microsoft.com/dotnet/sdk:10.0 AS build
WORKDIR /src

# Restore dependencies first (cached layer)
COPY src/Protein.Engine/Protein.Engine.csproj src/Protein.Engine/
RUN dotnet restore src/Protein.Engine/Protein.Engine.csproj

# Copy source and publish
COPY src/Protein.Engine/ src/Protein.Engine/
RUN dotnet publish src/Protein.Engine/Protein.Engine.csproj \
    -c Release \
    -o /app/publish \
    --no-restore

# ─────────────────────────────────────────────────────────────────────────────
# Stage 2: Runtime
# ─────────────────────────────────────────────────────────────────────────────
FROM mcr.microsoft.com/dotnet/runtime:10.0 AS runtime
WORKDIR /app

# Copy published binary
COPY --from=build /app/publish .

# Copy benchmark data and structure files
# (separate layer so code changes don't invalidate data cache)
COPY data/ /app/data/

# Non-secret defaults (URI and username are not credentials)
# NEO4J_PASSWORD must be supplied at runtime via docker-compose or -e flag
ENV NEO4J_URI=bolt://neo4j:7687
ENV NEO4J_USERNAME=neo4j

ENTRYPOINT ["dotnet", "Protein.Engine.dll"]
CMD ["--help"]
