using Infrastructure;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using System;
using System.Collections.Generic;
using System.Linq;
using Utility;
using static Org.BouncyCastle.Asn1.Cmp.Challenge;

namespace DataAccess
{
    public class DbInitializer : IDbInitializer
    {
        private readonly ApplicationDbContext _db;
        private readonly UserManager<AppUser> _userManager;
        private readonly RoleManager<IdentityRole> _roleManager;

        public DbInitializer(ApplicationDbContext db, UserManager<AppUser> userManager, RoleManager<IdentityRole> roleManager)
        {
            _db = db;
            _userManager = userManager;
            _roleManager = roleManager;
        }

        public void Initialize()
        {
            _db.Database.EnsureCreated();

            // Check if data has already been seeded
            if (_db.Rentals.Any())
            {
                return;
            }

            // Seeding amenities
            var amenities = new List<Amenity>
            {
                new Amenity { Name = "Hot Tub"},
                new Amenity { Name = "Big TV"},
                new Amenity { Name = "NES"},
                new Amenity { Name = "WiFi" },
                new Amenity { Name = "Swimming Pool" },
                new Amenity { Name = "Gym" },
                new Amenity { Name = "Parking" },
                new Amenity { Name = "Air Conditioning" },
                new Amenity { Name = "Pet Friendly" },
                new Amenity { Name = "Laundry Facilities" },
                new Amenity { Name = "Balcony" },
                new Amenity { Name = "Fireplace" }
            };

            foreach (var a in amenities)
            {
                _db.Amenities.Add(a);
            }
            _db.SaveChanges();

            // Seeding beds
            var beds = new List<Bed>
            {
                new Bed { BedType = "Tatami Mat"},
                new Bed { BedType = "Bunk Bed"},
                new Bed { BedType = "Twin"},
                new Bed { BedType = "Queen"},
                new Bed { BedType = "King"}
            };

            foreach (var m in beds)
            {
                _db.Beds.Add(m);
            }
            _db.SaveChanges();


            // Seeding roles and admin user
            _roleManager.CreateAsync(new IdentityRole(SD.AdminRole)).GetAwaiter().GetResult();
            _roleManager.CreateAsync(new IdentityRole(SD.CustomerRole)).GetAwaiter().GetResult();
            _roleManager.CreateAsync(new IdentityRole(SD.HostRole)).GetAwaiter().GetResult();



            _userManager.CreateAsync(new AppUser
            {
                UserName = "admin@bnb.com",
                Email = "admin@bnb.com",
                FName = "admin Fname",
                LName = "admin Lname",
                PhoneNumber = "1234567890",
                StreetAddress = "1234567890",
                State = "That one place",
                City = "That city in that one place",
                PostalCode = "12345",
            }, "Letmein1!").GetAwaiter().GetResult();

            AppUser user = _db.AppUser.FirstOrDefault(u => u.Email == "admin@bnb.com");
            _userManager.AddToRoleAsync(user, SD.AdminRole).GetAwaiter().GetResult();


            _userManager.CreateAsync(new AppUser
            {
                UserName = "host@bnb.com",
                Email = "host@bnb.com",
                FName = "host Fname",
                LName = "host Lname",
                PhoneNumber = "1234567890",
                StreetAddress = "1234567890",
                State = "That one place",
                City = "That city in that one place",
                PostalCode = "12345",
            }, "Letmein1!").GetAwaiter().GetResult();

            AppUser host = _db.AppUser.FirstOrDefault(u => u.Email == "host@bnb.com");
            _userManager.AddToRoleAsync(host, SD.HostRole).GetAwaiter().GetResult();


            // Seeding roles and Customer user
            _userManager.CreateAsync(new AppUser
            {
                UserName = "cust@bnb.com",
                Email = "cust@bnb.com",
                FName = "cust Fname",
                LName = "cust Lname",
                PhoneNumber = "0987654321",
                StreetAddress = "0987654321",
                State = "That other place",
                City = "That city in that other place",
                PostalCode = "54321",
            }, "Letmein1!").GetAwaiter().GetResult();

            AppUser cust = _db.AppUser.FirstOrDefault(u => u.Email == "cust@bnb.com");
            _userManager.AddToRoleAsync(cust, SD.CustomerRole).GetAwaiter().GetResult();



            // ------------Mass 100 User Seed--------------------------------------------------------------------------------------

            for (int i = 1; i <= 100; i++)
            {
                string email = $"user{i}@bnb.com";
                string role = i < 5 ? SD.HostRole : SD.CustomerRole;

                _userManager.CreateAsync(new AppUser
                {
                    UserName = email,
                    Email = email,
                    FName = $"User{i} Fname",
                    LName = $"User{i} Lname",
                    PhoneNumber = "1234567890",
                    StreetAddress = "1234567890",
                    State = "That one place",
                    City = "That city in that one place",
                    PostalCode = "12345",
                }, "Letmein1!").GetAwaiter().GetResult();

                AppUser newUser = _db.AppUser.FirstOrDefault(u => u.Email == email);
                _userManager.AddToRoleAsync(newUser, role).GetAwaiter().GetResult();
            }

            //---------------------------------------------------------------------------------------------------------------------





            var owner = _db.AppUser.FirstOrDefault(u => u.Email == "host@bnb.com");

            // Seeding rentals with associated rental amenities, beds, prices, images and review
            var rentals = new List<Rental>
            {

                new Rental {
                    OwnerName = owner,
                    Name = "Cozy Apartment",
                    Description = "A lovely apartment perfect for couples",
                    Address = "123 Main Street",
                    City = "New York",
                    State = "NY",
                    Zip = "10001",
                    RentalType = "Apartment",
                    MaxOcupants = 2,
                    Bathrooms = 1,
                    Beds = 3,
                    Published = true,
                    RentalAmenities = new List<RentalAmenity>
                    {
                        new RentalAmenity { Amenity = amenities[0] }, // Hot Tub
                        new RentalAmenity { Amenity = amenities[3] }  // WiFi
                    },
                    RentalBeds = new List<RentalBed>
                    {
                        new RentalBed { Bed = beds[0] }, // Tatami Mat
                        new RentalBed { Bed = beds[2] }  // Twin
                    },
                    Prices = new List<Price>
                    {
                        new Price
                        {
                            PricePerNight = 100.00,
                            PriceDateStart = DateTime.Now,
                            PriceDateEnd = DateTime.Now.AddDays(30)
                        }
                    },
                    Images = new List<Image> 
                    {
                        new Image
                        {
                            RentalImage = "/images/rentalimg/cozy-apartment.jpg",
                            PrimaryImage = true
                        },
                        new Image
                        {
                            RentalImage = "/images/rentalimg/cozy-apartment2.jpg",
                            PrimaryImage = false
                        },
                        new Image
                        {
                            RentalImage = "/images/rentalimg/cozy-apartment3.jpg",
                            PrimaryImage = false
                        },
                    },
                    Reviews = new List<Review> 
                    {
                        new Review
                        {
                            Reviewer = cust, 
                            body = "This apartment was amazing! Highly recommended."
                        }
                    }
                },
                new Rental {
                    OwnerName = owner,
                    Name = "Spacious House",
                    Description = "A spacious house with a backyard",
                    Address = "456 Elm Avenue",
                    City = "Los Angeles",
                    State = "CA",
                    Zip = "90001",
                    RentalType = "House",
                    MaxOcupants = 6,
                    Bathrooms = 2,
                    Beds = 6,
                    Published = true,
                    RentalAmenities = new List<RentalAmenity>
                    {
                        new RentalAmenity { Amenity = amenities[1] }  // Big TV
                    },
                    RentalBeds = new List<RentalBed>
                    {
                        new RentalBed { Bed = beds[1] },
                        new RentalBed { Bed = beds[2] },
                        new RentalBed { Bed = beds[2] },
                        new RentalBed { Bed = beds[3] },
                        new RentalBed { Bed = beds[3] },
                        new RentalBed { Bed = beds[4] },
                    },
                    Prices = new List<Price>
                    {
                        new Price
                        {
                            PricePerNight = 75.00,
                            PriceDateStart = DateTime.Now,
                            PriceDateEnd = DateTime.Now.AddDays(25)
                        }
                    },
                    Images = new List<Image>
                    {
                        new Image
                        {
                            RentalImage = "/images/rentalimg/Spacious House.jpg",
                            PrimaryImage = true
                        },
                        new Image
                        {
                            RentalImage = "/images/rentalimg/Spacious House2.jpg",
                            PrimaryImage = false
                        },
                        new Image
                        {
                            RentalImage = "/images/rentalimg/Spacious House3.jpg",
                            PrimaryImage = false
                        },
                    },
                    Reviews = new List<Review>
                    {
                        new Review
                        {
                            Reviewer = cust,
                            body = "A super cool review for this house"
                        }
                    }
                },
                new Rental {
                    OwnerName = owner,
                    Name = "Luxury Condo",
                    Description = "A luxurious condo with a great view",
                    Address = "789 Oak Boulevard",
                    City = "Miami",
                    State = "FL",
                    Zip = "33101",
                    RentalType = "Condo",
                    MaxOcupants = 4,
                    Bathrooms = 2,
                    Beds = 3,
                    Published = true,
                    RentalAmenities = new List<RentalAmenity>
                    {
                        new RentalAmenity { Amenity = amenities[0] },
                        new RentalAmenity { Amenity = amenities[1] },
                        new RentalAmenity { Amenity = amenities[3] }
                    },
                    RentalBeds = new List<RentalBed>
                    {
                        new RentalBed { Bed = beds[0] },
                        new RentalBed { Bed = beds[1] },
                        new RentalBed { Bed = beds[4] }
                    },
                    Prices = new List<Price>
                    {
                        new Price
                        {
                            PricePerNight = 2000.00,
                            PriceDateStart = DateTime.Now,
                            PriceDateEnd = DateTime.Now.AddDays(365)
                        }
                    },
                    Images = new List<Image>
                    {
                        new Image
                        {
                            RentalImage = "/images/rentalimg/LuxaryCondo.jpg",
                            PrimaryImage = true
                        },
                        new Image
                        {
                            RentalImage = "/images/rentals/IMAGEHERE.jpg",
                            PrimaryImage = false
                        },
                    },
                    Reviews = new List<Review>
                    {
                        new Review
                        {
                            Reviewer = cust,
                            body = "A super cool review for this Condo"
                        }
                    }
                }
            };



            //---------------------------------------------------------------------------------------------------------------- Add 100  more random rentals
            Random random = new Random();

            // Fetch all hosts from the database based on the HostRole and add them to newHosts
            var hostIds = _userManager.GetUsersInRoleAsync(SD.HostRole).GetAwaiter().GetResult().Select(u => u.Id).ToList();
            var newHosts = _db.AppUser.AsEnumerable().Where(u => hostIds.Contains(u.Id)).ToList();

            //Mass Rental Seed
            for (int i = 0; i < 100; i++)
            {
                var rowner = newHosts[random.Next(newHosts.Count)];

                var extraRental = new Rental
                {
                    OwnerName = rowner,
                    Name = $"Random Rental {i}",
                    Description = $"Description for Random Rental {i}",
                    Address = $"Random Address {i}",
                    City = "Random City",
                    State = "Random State",
                    Zip = $"Random Zip {i}",
                    RentalType = "Random Rental Type", // Adjust this based on your types
                    MaxOcupants = random.Next(1, 10),
                    Bathrooms = random.Next(1, 4),
                    Beds = random.Next(1, 6),
                    Published = true,
                    RentalAmenities = new List<RentalAmenity>
                    {
                        new RentalAmenity { Amenity = amenities[random.Next(amenities.Count)] },
                        new RentalAmenity { Amenity = amenities[random.Next(amenities.Count)] }
                    },
                    RentalBeds = new List<RentalBed>
                    {
                        new RentalBed { Bed = beds[random.Next(beds.Count)] },
                        new RentalBed { Bed = beds[random.Next(beds.Count)] }
                    },
                    Prices = new List<Price>
                    {
                        new Price
                        {
                            PricePerNight = random.Next(50, 500),
                            PriceDateStart = DateTime.Now,
                            PriceDateEnd = DateTime.Now.AddDays(random.Next(10, 90))
                        }
                    },
                    Images = new List<Image>
                    {
                        new Image
                        {
                            RentalImage = $"images/rentalimg/house{i+1}.jpg",
                            PrimaryImage = true
                        },
                        new Image
                        {
                            RentalImage = $"images/rentalimg/house{i+1}.jpg",
                            PrimaryImage = false
                        }
                    },
                    Reviews = new List<Review>
                    {
                        new Review
                        {
                            Reviewer = cust,
                            body = $"A review for rental{i}"
                        }
                    }
                };

                rentals.Add(extraRental);
            }

            //----------------------------------------------------------------------------------------------------------------End 100 more random rentals

            foreach (var r in rentals)
            {
                _db.Rentals.Add(r);
            }
            _db.SaveChanges();







            //Test seed of reservation 
            var cozyApartment = _db.Rentals.FirstOrDefault(r => r.Name == "Cozy Apartment");

            if (cozyApartment != null)
            {
                // Create a reservation for the Cozy Apartment
                var reservation = new Reservation
                {
                    Rental = cozyApartment.Id,
                    IndividualTotal = 100.00, 
                    SalesTax = 10.00, 
                    GrandTotal = 110.00,
                    Guests = 2, 
                    Status = 2, // something like 1 = Pending, 2 = Approved, 3 = Cancelled
                    BookedStart = DateTime.Now.AddDays(5), 
                    BookedEnd = DateTime.Now.AddDays(10),
                    AppUser = cust
                };

                _db.Reservation.Add(reservation);
                _db.SaveChanges();
            }
        }
    }
}