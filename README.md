# Solar-System-Graphics
Computer Graphics project

This project uses multiple graphics libraries and tools such as GLUT, GLEW and OpenGL. It also processes XML files using the library TinyXML to parse these file types.

The purpose of this project is to learn how to compute/draw primitive and complex graphics models, such as a teapot and another complex geometric figures.

Two core components define this project: generator and motor. Allowing to draw custom environments by using the multiple models attached with textures and animations like translations and movements along predefined paths.

The generator is responsible for generate files that defines a geometric figure. On these files are printed all the vertices necessary to the motor draw the figure.

The motor parses the generated model files and compute its graphic modeling by joining edges between all defined vertices.

Moreover, bezier patches, catmull-Rom curves, textures appliance, movements on figures and XML parsing and manipulation is also explored on this project.
