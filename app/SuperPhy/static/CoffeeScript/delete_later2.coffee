###

Delete Later Module 2
-------------------

This is just an example of running browserify on a project to:
1. Compile coffeescript
2. Use of require and exports
3. Bundle and generate single .js
4. Add source maps for debugging



author: Matt Whtieside
copyright: Copyright 2016, Public Health Agency of Canada
license: ASL
version: 2.0
maintainer = Matt Whiteside
email: matthew.whiteside at phac-aspc.gc.ca

###

# Function that is exported from this module
module.exports = (name) ->
  return 'Hello ' + name + '!'

