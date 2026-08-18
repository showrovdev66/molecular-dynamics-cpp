# GitHub Setup Notes

## Suggested repository name

```text
molecular-dynamics-cpp
```

## Suggested GitHub About description

```text
Three learning projects in C++: 2D Lennard-Jones MD, cell-list force optimization, and 3D reverse-NEMD thermal conductivity.
```

## Suggested topics

```text
cplusplus
molecular-dynamics
computational-physics
lennard-jones
neighbor-list
cell-list
nemd
thermal-conductivity
simulation
```

## What should be uploaded

Upload the complete contents of this prepared folder.

Do not upload old Windows `.exe` files or temporary context/instruction notes. They are intentionally excluded from this repository.

## First Git upload

Open a terminal inside the repository folder and run:

```bash
git init
git add .
git commit -m "Add molecular dynamics learning projects"
git branch -M main
git remote add origin YOUR_GITHUB_REPOSITORY_URL
git push -u origin main
```

Replace `YOUR_GITHUB_REPOSITORY_URL` with the URL GitHub gives you after you create the empty repository.

## Recommended visibility

Use **Public** if you want to show the projects as part of a portfolio.

## Before publishing

Check these points:

- the root `README.md` appears correctly on GitHub;
- the image inside Project 1 loads;
- no `.exe` file is committed;
- no temporary context or instruction file is committed;
- all three source files compile locally;
- Project 1's Python plot script runs;
- benchmark timing is described as machine-dependent;
- the NEMD project is described as an educational reduced-size implementation, not an exact reproduction of the publication.
