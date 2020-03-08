from fluid_cube import FluidCube

if __name__ == '__main__':
    cube = FluidCube(100, 1, 1, 0.01)
    cube.add_density(50, 50, 50, 100)
    cube.add_velocity(50, 50, 50, 1, -1, 1)
    cube.fluid_step()
    print(cube.density)
