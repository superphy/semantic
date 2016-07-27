/*jslint node: true */

'use strict';
var gulp = require('gulp');
var coffee = require('gulp-coffee');
var gutil = require('gutil');
var uglify = require('gulp-uglify')
var rename = require('gulp-rename');
var less = require('gulp-less');
var minifyCSS = require('gulp-minify-css');
var concat = require('gulp-concat');
var flatten = require('gulp-flatten');
var order = require('gulp-order');
var print = require('gulp-print');
var browserify  = require('browserify');
var source      = require('vinyl-source-stream'); //to 'rename' your resulting file
var buffer      = require('vinyl-buffer'); // to transform the browser
var sourcemaps  = require('gulp-sourcemaps');
var coffeeify  = require('coffeeify');

var ordering = [
  "__init__.coffee",
  "Pages/Page.coffee",
  "Components/**/*.coffee",
  "Pages/**/*.coffee",
  "Models/**/*.coffee",
  "!route.coffee",
  "**/*.coffee",
  "route.coffee"
];

var config = {
    coffee: {
        src: 'CoffeeScript/delete_later.coffee',       // Entry point
        outputDir: './js/',  // Directory to save bundle to
        mapDir: './maps/',      // Subdirectory to save maps to
        outputFile: 'superphy_bundle.js' // Name to use for bundle
    },
};

// This method makes it easy to use common bundling options in different tasks
function bundle (bundler) {

    // Add options to add to "base" bundler passed as parameter
    bundler
      .bundle()                                                        // Start bundle
      .pipe(source(config.coffee.src))                                 // Entry point
      .pipe(buffer())                                                  // Convert to gulp pipeline
      .pipe(rename(config.coffee.outputFile))                          // Rename output
      .pipe(sourcemaps.init({ loadMaps : true, debug: true }))  // Strip inline source maps
      .pipe(uglify())  // Minify
      .pipe(sourcemaps.write(config.coffee.mapDir))    // Save source maps to own directory
      .pipe(gulp.dest(config.coffee.outputDir))        // Save 'bundle'
      //.pipe(livereload());                                       // Reload browser if relevant
}

gulp.task('bundle', function () {
    var bundler = browserify(config.coffee.src)  // Pass browserify the entry point
                    .transform(coffeeify)      //  Chain transformations: First, coffeeify . . .
                              

    bundle(bundler);  // Chain other options -- sourcemaps, rename, etc.
})

gulp.task('coffee_to_js',  function() {
	gulp
    .src("CoffeeScript/**/*.coffee")
    .pipe(order(ordering))
    .pipe(print(function(filepath) {
      return "build: " + filepath;
    }))
  	.pipe(flatten())
  	.pipe(concat('superphy.coffee'))
  	.pipe(coffee())
    //.pipe(uglify())
    .pipe( rename( { extname: '.js'} ) )
    .pipe(gulp.dest('./js'));
});

gulp.task('less_to_css', function() {
  return gulp.src(['./less/*.less'])
    .pipe(less())
    .pipe(gulp.dest('./css/'));
});

gulp.task('minify_css', function() {
  gulp.src('./css/*.css')
    .pipe(concat('all.css'))
    .pipe(minifyCSS())
    .pipe(rename({
      extname: '.min.css'
    }))
    .pipe(gulp.dest('./css/'))
});

gulp.task('default', ['less_to_css'], function() {
  gulp.start('coffee_to_js', 'minify_css');
});