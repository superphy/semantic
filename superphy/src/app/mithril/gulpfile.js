/*jslint node: true */

'use strict';
var gulp = require('gulp');
var coffee = require('gulp-coffee');
var gutil = require('gutil');
var uglify = require('gulp-uglify');
var rename = require('gulp-rename');
var less = require('gulp-less');
var minifyCSS = require('gulp-minify-css');
var concat = require('gulp-concat');
var flatten = require('gulp-flatten');
var streamqueue = require('streamqueue');

gulp.task('coffee_to_js',  function() {
	return streamqueue ({ objectMode: true},
		//Templates are inherited by both Components and Pages
    gulp.src(['./coffee/Templates.coffee']),
    //Components are used in Pages
    gulp.src(['./coffee/Components/*.coffee']),
    //Pages are used in routes. Routes are found in root folder
    gulp.src(['./coffee/Pages/*.coffee']),
    gulp.src(['./coffee/Pages/*/*.coffee']),
    //All other subfolders are imported before the root
    gulp.src(['./coffee/*/*.coffee']),
    //All other root files are imported before __init__
    gulp.src(['./coffee/!(__init__)*.coffee']),
    //By including __init__ last, every other class is before this in the single page aplication.
    gulp.src(['./coffee/__init__.coffee'])
	)
  	   	.pipe(flatten())
		.pipe(concat('all.coffee'))
		.pipe(coffee())
      	//.pipe(uglify())
        .pipe(rename({
        	extname: '.min.js'
        }))
        .pipe(gulp.dest('./js')
    );
});

gulp.task('less_to_css', function() {
  return gulp.src(['./less/*.less'])
      .pipe(less())
      .pipe(gulp.dest('./css/'));
});

gulp.task('minify_css', function() {
  return gulp.src('./css/*.css')
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