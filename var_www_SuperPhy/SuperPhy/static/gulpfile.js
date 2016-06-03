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