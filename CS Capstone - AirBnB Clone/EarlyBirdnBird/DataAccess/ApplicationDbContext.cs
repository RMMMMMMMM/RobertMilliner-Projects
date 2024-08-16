using Microsoft.EntityFrameworkCore;
using Infrastructure.Models;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Microsoft.AspNetCore.Identity.EntityFrameworkCore;

namespace DataAccess
{
    public class ApplicationDbContext : IdentityDbContext
    {
        public ApplicationDbContext(DbContextOptions<ApplicationDbContext> options) : base(options)
        {
        }

        public DbSet<Amenity> Amenities { get; set; }

        public DbSet<Bed> Beds { get; set; }

        public DbSet<Block> Blocked { get; set; }


        public DbSet<Image> Images { get; set; }

        public DbSet<Price> Prices { get; set; }

        public DbSet<Rental> Rentals { get; set; }

        public DbSet<Reservation> Reservation { get; set; }

        public DbSet<AppUser> AppUser { get; set; }

        public DbSet<RentalAmenity> RentalAmenities { get; set; }

        public DbSet<UserMail> UserMails { get; set; }

        public DbSet<Review> Reviews { get; set; }


    }
}
